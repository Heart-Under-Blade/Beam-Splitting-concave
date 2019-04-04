#include "Tracks.h"

#include <iostream>

#include "Beam.h"

Tracks::Tracks()
{
	tree = new TrackNode(-1);
}

int Tracks::FindGroupByTrackId(const IdType &trackId) const
{
	for (int i = 0; i < size(); ++i)
	{
		for (int j = 0; j < (*this)[i].size; ++j)
		{
			if ((*this)[i].arr[j] == trackId)
			{
				return (*this)[i].groupID;
			}
		}
	}

	if (size() == 0)
	{
		return 0;
	}

	return -1;
}

void Tracks::FillTrackTree(const std::vector<std::vector<int>> &tracks)
{
	auto currNode = tree;

	for (auto &t : tracks)
	{
		for (int i = 0; i < t.size(); ++i)
		{
			int id = t[i];
			auto node = currNode->FindNode(id);

			if (node != nullptr)
			{
				currNode = node;
			}
			else
			{
				currNode = currNode->AddChild(id);
			}

			if (i == t.size()-1)
			{
				currNode->isCompleted = true;
			}
		}

		currNode = tree;
	}
}

void Tracks::CreateGroupsForUngroupedTracks(const TrackGroup &buffGroup)
{
	if (buffGroup.size != 0) // добавляем треки без группы в отдельные группы
	{
		for (int i = 0; i < buffGroup.size; ++i)
		{
			TrackGroup newGroup;
			newGroup.arr[newGroup.size++] = buffGroup.arr[i];
			newGroup.groupID = size();
			push_back(newGroup);
		}
	}
}

int Tracks::ImportTrack(char *buff, vector<int> &track)
{
	int groupIndex = -1;
	char *ptr, *trash;
	ptr = strtok(buff, " ");

	while (ptr != NULL)
	{
		if (ptr[0] == ':')
		{
			ptr = strtok(NULL, " ");
			groupIndex = strtol(ptr, &ptr, 10);

			if (groupIndex >= size())
			{
				for (int i = size(); i <= groupIndex; ++i)
				{
					(*this).push_back(TrackGroup());
				}
			}

			(*this)[groupIndex].groupID = groupIndex;
			break;
		}

		int tmp = strtol(ptr, &trash, 10);
		track.push_back(tmp);
		ptr = strtok(NULL, " ");
	}

	return groupIndex;
}

IdType Tracks::ComputeTrackId(const vector<int> &track, int nFacets)
{
	IdType trackId = 0;

	for (int t : track)
	{
		trackId += (t + 1);
		trackId *= (nFacets + 1);
	}

	return trackId;
}

void Tracks::ImportTracks(int nFacets, const std::string &filename)
{
	std::vector<std::vector<int>> tracks;

	m_nFacets = nFacets;
	const int bufSize = 1024;
	std::ifstream trackFile(filename, std::ios::in);

	if (!trackFile.is_open())
	{
		std::cerr << "Track file not found" << std::endl;
		throw std::exception();
	}

	char *buff = (char*)malloc(sizeof(char) * bufSize);

	TrackGroup buffGroup;

	while (!trackFile.eof())
	{
		trackFile.getline(buff, bufSize);

		vector<int> track;
		int groupIndex = ImportTrack(buff, track);

		auto trackId = ComputeTrackId(track, nFacets);

		if (groupIndex >= 0)
		{
			(*this)[groupIndex].arr[(*this)[groupIndex].size++] = trackId;
		}
		else
		{
			buffGroup.arr[buffGroup.size++] = trackId;
		}

		tracks.push_back(track);

		track.clear();
	}

	trackFile.close();

	CreateGroupsForUngroupedTracks(buffGroup);
	FillTrackTree(tracks);
}

void Tracks::RecoverTrack(const Beam &beam, std::vector<int> &track)
{
	RecoverTrack(m_nFacets, beam, track);
}

void Tracks::RecoverTrack(int nFacets, const Beam &beam, std::vector<int> &track)
{
	int coef = nFacets + 1;
	std::vector<int> tmp_track;

	auto tmpId = beam.id/coef;

	for (int i = 0; i <= beam.actNo; ++i)
	{
#ifdef _DEBUG // DEB
		int tmp = (tmpId%coef);
#else
		int tmp = (tmpId%coef).toInt();
#endif
		tmpId -= tmp;
		tmpId /= coef;
		tmp -= 1;
		tmp_track.push_back(tmp);
	}

	for (int i = tmp_track.size()-1; i >= 0; --i)
	{
		track.push_back(tmp_track.at(i));
	}
}

std::string Tracks::TrackToStr(const std::vector<int> &track)
{
	std::string str;

	for (int p : track)
	{
		str += std::to_string(p) + ' ';
	}

	return str;
}
