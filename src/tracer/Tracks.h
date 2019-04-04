#pragma once

#include <cstring>
#include <string>
#include <vector>

#include "BigInteger.hh"
#include "TrackTree.h"

#define MAX_GROUP_NUM	1024

#ifdef _DEBUG // DEB
typedef long long IdType;
#else
typedef BigInteger IdType;
#endif

class Beam;
/// REF OPT TODO: сохранять треки вместе в id в виде массива интов,
/// чтобы не конвертировать, а просто искать их по id
class TrackGroup
{
public:
	int groupID;
	IdType arr[MAX_GROUP_NUM];
	int size = 0;

	std::string CreateGroupName() const
	{
		std::string subname;
		subname += "gr_" + std::to_string(groupID);
		return subname;
	}
};

class Tracks : public std::vector<TrackGroup>
{
public:
	Tracks();
	Tracks(int nFacets) {m_nFacets = nFacets;}
	int FindGroupByTrackId(const IdType &trackId) const;

	void ImportTracks(int nFacets, const std::string &filename);
	void RecoverTrack(const Beam &beam, std::vector<int> &track);
	static void RecoverTrack(int nFacets, const Beam &beam, std::vector<int> &track);
	static std::string TrackToStr(const std::vector<int> &track);

	bool shouldComputeTracksOnly;
	TrackNode *tree;

	IdType ComputeTrackId(const vector<int> &track, int nFacets);

private:
	int m_nFacets;

	void FillTrackTree(const std::vector<std::vector<int>> &tracks);
	void CreateGroupsForUngroupedTracks(const TrackGroup &buffGroup);
	int ImportTrack(char *buff, vector<int> &track);
};
