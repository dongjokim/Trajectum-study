#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include "AliJBaseEventHeader.h"
#include "AliJBaseTrack.h"

static const double kPerc[] = {0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
static const size_t kNPerc = sizeof(kPerc) / sizeof(kPerc[0]);

static std::vector<double> LoadCentRow(const std::string &csv, int paramID) {
	std::ifstream file(csv);
	if (!file)
		return {};
	for (int i = 0; i < paramID; ++i)
		file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::string line;
	if (!std::getline(file, line))
		return {};
	std::vector<double> bins;
	std::stringstream ss(line);
	std::string tok;
	while (std::getline(ss, tok, ','))
		bins.push_back(std::stod(tok));
	return bins;
}

// Percentile table from this sample: high Nch = central (0%).
static std::vector<double> EdgesFromMult(std::vector<int> nch) {
	if (nch.empty())
		return {};
	std::sort(nch.begin(), nch.end(), std::greater<int>());
	std::vector<double> edges;
	edges.reserve(kNPerc);
	const double n = static_cast<double>(nch.size() - 1);
	for (size_t i = 0; i < kNPerc; ++i) {
		const double idx = std::min(n, std::max(0.0, kPerc[i] / 100.0 * n));
		const size_t lo = static_cast<size_t>(idx);
		const size_t hi = std::min(nch.size() - 1, lo + 1);
		const double frac = idx - static_cast<double>(lo);
		edges.push_back(nch[lo] * (1.0 - frac) + nch[hi] * frac);
	}
	return edges;
}

static void WriteCentCsv(const std::string &path, const std::vector<double> &edges) {
	std::ofstream out(path);
	out << std::fixed << std::setprecision(4);
	for (size_t i = 0; i < edges.size(); ++i) {
		if (i)
			out << ',';
		out << edges[i];
	}
	out << '\n';
}

static float CentralityFromMult(double dnch, const std::vector<double> &edges) {
	if (edges.size() < 2)
		return 0.f;
	const size_t n = std::min(edges.size(), kNPerc);
	if (dnch >= edges[0])
		return 0.f;
	for (size_t i = 0; i + 1 < n; ++i) {
		if (dnch < edges[i + 1])
			continue;
		const double frac = (edges[i] - dnch) / (edges[i] - edges[i + 1] + 1e-12);
		return static_cast<float>(std::max(kPerc[i] + frac * (kPerc[i + 1] - kPerc[i]), 0.0));
	}
	return 99.f;
}

struct OscarLine {
	double px, py, pz, p0;
	long long pdg, id;
	int charge;
};

static bool ParseParticle(const std::string &line, OscarLine &p) {
	std::stringstream ss(line);
	double t, x, y, z, mass;
	if (!(ss >> t >> x >> y >> z >> mass >> p.p0 >> p.px >> p.py >> p.pz >> p.pdg >> p.id >> p.charge))
		return false;
	return p.charge != 0;
}

static double Eta(double px, double py, double pz) {
	const double p = std::sqrt(px * px + py * py + pz * pz);
	if (p <= 0.0)
		return 0.0;
	return 0.5 * std::log((p + pz) / (p - pz + 1e-12));
}

template <typename Fn>
static int ForEachOscarEvent(const std::string &path, Fn &&fn) {
	std::ifstream in(path);
	if (!in)
		return -1;
	std::vector<OscarLine> ev;
	int eventID = -1;
	bool inEvent = false;
	int nEv = 0;
	std::string line;
	auto flush = [&]() {
		if (!inEvent)
			return;
		fn(eventID, ev);
		++nEv;
		ev.clear();
		inEvent = false;
	};
	while (std::getline(in, line)) {
		if (line.empty())
			continue;
		if (line[0] == '#') {
			if (line.rfind("# event", 0) == 0) {
				const bool isEnd = line.find("end") != std::string::npos;
				if (isEnd) {
					flush();
				} else {
					flush();
					++eventID;
					inEvent = true;
				}
			}
			continue;
		}
		if (!inEvent) {
			inEvent = true;
			++eventID;
		}
		OscarLine p;
		if (ParseParticle(line, p))
			ev.push_back(p);
	}
	if (inEvent)
		flush();
	return nEv;
}

static int Nch05(const std::vector<OscarLine> &ev) {
	int n = 0;
	for (const auto &p : ev) {
		if (std::fabs(Eta(p.px, p.py, p.pz)) <= 0.5)
			++n;
	}
	return n;
}

int main(int argc, char **argv) {
	if (argc < 3) {
		std::cerr << "Usage: " << argv[0]
			  << " <oscar> <out.root> [cent.csv] [paramID]\n"
			  << "  With no cent.csv, percentiles are measured from this Oscar sample.\n";
		return 1;
	}
	const std::string inPath = argv[1];
	const std::string outPath = argv[2];
	std::vector<double> centEdges;
	if (argc >= 4)
		centEdges = LoadCentRow(argv[3], argc >= 5 ? std::atoi(argv[4]) : 0);

	std::vector<int> nchAll;
	if (centEdges.size() < 2) {
		if (ForEachOscarEvent(inPath, [&](int, const std::vector<OscarLine> &ev) {
			    nchAll.push_back(Nch05(ev));
		    }) < 0) {
			std::cerr << "Cannot open " << inPath << "\n";
			return 1;
		}
		centEdges = EdgesFromMult(nchAll);
		const std::string csvPath = outPath + ".cent.csv";
		WriteCentCsv(csvPath, centEdges);
		std::cout << "Centrality from this sample (" << nchAll.size()
			  << " events), Nch(|eta|<=0.5) at percentiles\n";
		std::cout << "  0,1,5,10,20,30,40,50,60,70,80,90,100% :";
		for (double e : centEdges)
			std::cout << ' ' << e;
		std::cout << "\n  wrote " << csvPath << "\n";
	}

	TFile out(outPath.c_str(), "RECREATE");
	TTree jTree("jTree", "JCORRAN tree from Trajectum/SMASH Oscar");
	auto *tracks = new TClonesArray("AliJBaseTrack", 4000);
	auto *headers = new TClonesArray("AliJBaseEventHeader", 1);
	jTree.Branch("JTrackList", &tracks);
	jTree.Branch("JEventHeaderList", &headers);

	int written = 0;
	ForEachOscarEvent(inPath, [&](int eventID, const std::vector<OscarLine> &ev) {
		tracks->Clear();
		for (const auto &p : ev) {
			const int idx = tracks->GetEntriesFast();
			new ((*tracks)[idx]) AliJBaseTrack(p.px, p.py, p.pz, p.p0, static_cast<int>(p.id),
							   static_cast<Short_t>(p.pdg),
							   static_cast<Char_t>(p.charge));
		}
		const int nch05 = Nch05(ev);
		const float cent = CentralityFromMult(nch05, centEdges);
		headers->Clear();
		auto *hdr = new ((*headers)[0]) AliJBaseEventHeader(eventID, cent, 0.f);
		hdr->SetXVertexMC(static_cast<float>(nch05));
		hdr->SetYVertexMC(static_cast<float>(nch05));
		jTree.Fill();
		++written;
	});

	jTree.Write();
	out.Close();
	std::cout << "Wrote " << written << " events to " << outPath << "\n";
	return written > 0 ? 0 : 1;
}
