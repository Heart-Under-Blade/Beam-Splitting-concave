#pragma once

#include "argparse.hpp"

class ArgParser : public ArgumentParser
{
public:
	template <class T>
	T getArgValue(const std::string &arg)
	{
		std::string arg_str = retrieve<std::string>(arg);
		return argToValue<T>(arg_str);
	}

	template <class T>
	T argToValue(const std::string &arg) const
	{
		bool ok = true;
		char *end;
		T val;

		if (typeid(T) == typeid(int))
		{
			val = (T)strtol(arg.c_str(), &end, 10);
		}
		else if (typeid(T) == typeid(double))
		{
			val = (T)strtod(arg.c_str(), &end);
		}
		else
		{
			ok = false;
		}

		if (strlen(end) != 0 || !ok)
		{
			throw std::string("Some argument is incorrect.");
		}

		return val;
	}

};
