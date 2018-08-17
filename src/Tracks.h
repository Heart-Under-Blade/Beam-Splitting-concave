#pragma once

#include <cstring>
#include <string>
#include <vector>

#include "BigInteger.hh"

#define MAX_GROUP_NUM	4096

class Beam;
/// REF OPT TODO: сохранять треки вместе в id в виде массива интов,
/// чтобы не конвертировать, а просто искать их по id
class TrackGroup
{
public:
	int groupID;
	BigInteger arr[MAX_GROUP_NUM];
	int size = 0;
	std::vector<std::vector<int>> tracks;

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
	int FindGroupByTrackId(const BigInteger &trackId) const;

	void ImportTracks(int nFacets, const std::string &filename);
	static void RecoverTrack(const Beam &beam, int facetNum,
							 std::vector<int> &track);

	bool shouldComputeTracksOnly;
};
