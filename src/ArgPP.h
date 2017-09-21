#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <typeinfo>
#include <cstring>

class ArgPP
{
public: // methods

	ArgPP() {}

	// TODO: write 'throw' in method's declarations
	void AddRule(const std::string &key, char valueNum = 0,
				 bool isOptional = false, const std::string &dependsOn = "")
	{
		Rule rule(valueNum, isOptional, dependsOn);
		AddArgRule(key, rule);
	}

	void AddRule(const std::string &key, const std::string &longKey, char valueNum = 0,
				 bool isOptional = false, const std::string &dependsOn = "")
	{
		if (longKey.size() > 1)
		{
			Rule rule1(valueNum, isOptional, dependsOn, longKey);
			AddArgRule(key, rule1);

			Rule rule2(valueNum, isOptional, dependsOn, key);
			AddArgRule(longKey, rule2);
		}
		else
		{
			std::string msg = "long key must consists of more than 1 chatacter";
			Error(msg);
		}
	}

	void Parse(int argc, const char *argv[])
	{
		m_args.clear();
		m_programName = argv[0];
		std::vector<std::string> rawArgs(argv+1, argv+argc);
		FillArgs(rawArgs);
		CheckRequiredArgs();
	}

	int GetIntValue(const std::string &key, size_t valueIndex = 0) const
	{
		std::string rawValue = FindArgValue(key, valueIndex);
		return ConvertTo<int>(rawValue);
	}

	double GetDoubleValue(const std::string &key, size_t valueIndex = 0) const
	{
		std::string rawValue = FindArgValue(key, valueIndex);
		return ConvertTo<double>(rawValue);
	}

	std::string GetStringValue(const std::string &key, size_t valueIndex = 0) const
	{
		return FindArgValue(key, valueIndex);
	}

	// TODO: несколько раз вызываешь - возвращает каждый раз сдед. значение аргумента
	int GetNextIntValue(const std::string &key) const
	{
		// ...
	}

	const std::string &GetProgramName() const
	{
		return m_programName;
	}

	bool Catched(const std::string &key) const
	{
		auto it = m_args.find(key);
		return it != m_args.end();
	}

	void Reset()
	{
		m_rules.clear();
		m_args.clear();
		m_programName.clear();
	}

private: // types

	struct Rule
	{
		Rule(char valueNum, bool isOptional, const std::string &dependsOn,
			 const std::string &secondaryBeam = "")
			: valueNum(valueNum),
			  isOptional(isOptional),
			  dependsOn(dependsOn),
			  secondaryKey(secondaryBeam)
		{
		}

		char valueNum;
		bool isOptional;
		std::string dependsOn;
		std::string secondaryKey;
	};

	typedef std::pair<std::string, Rule> NamedRule;
	typedef std::vector<std::string> Arg;
	typedef std::pair<std::string, Arg> NamedArg;

private: // fields

	std::map<std::string, Rule> m_rules;
	std::map<std::string, Arg> m_args;
	std::string m_programName;

private: // methods

	std::string FindArgValue(const std::string &key, size_t valueIndex) const
	{
		Arg arg = FindArg(key);

		if (valueIndex >= arg.size())
		{
			std::string msg = "cannot find value " + std::to_string(valueIndex)
					+ " of argument '" + arg[valueIndex] + '\'';
			Error(msg);
		}

		return arg[valueIndex];
	}

	void FillArgs(const std::vector<std::string> &rawArgs)
	{
		std::string key;
		size_t valueNum = 0;

		for (size_t i = 0; i < rawArgs.size(); ++i)
		{
			const std::string &rawArg = rawArgs[i];

			if (valueNum == 0) // read key
			{
				key = RetrieveKey(rawArg);
				Rule rule = FindRule(key);
				CheckDependency(key, rule);
				m_args.insert(NamedArg(key, Arg()));
				valueNum = rule.valueNum;
			}
			else // read key value (one of these)
			{
				FillArg(key, rawArg, i, valueNum);
			}
		}

		if (valueNum != 0)
		{
			std::string msg = "not enough values for argument with key '"
					+ key + "', must be " + std::to_string(valueNum) + " more";
			Error(msg);
		}
	}

	void FillArg(const std::string &key, const std::string &rawArg,
				 size_t &i, size_t &valueNum)
	{
		if (valueNum == '+') // more than one args
		{
			if (IsKey(rawArg[0]))
			{
				--i;
				valueNum = 0;

				if (!m_args[key].empty())
				{
					return; // TODO: kostyl'
				}
			}

			AddArgValue(key, rawArg);
		}
		else if (valueNum == '*') // one or more args
		{
			if (IsKey(rawArg[0]))
			{
				--i;
				valueNum = 0;
			}
			else
			{
				AddArgValue(key, rawArg);
			}
		}
		else // fixed number of args
		{
			AddArgValue(key, rawArg);
			--valueNum;
		}
	}

	const Rule &FindRule(const std::string &key)
	{
		auto it = m_rules.find(key);

		if (it == m_rules.end())
		{
			std::string msg = "rule for the argument with key '" + key
					+ "' not found";
			Error(msg);
		}

		return it->second;
	}

	const Arg &FindArg(const std::string &key) const
	{
		auto it = m_args.find(key);

		if (it == m_args.end())
		{
			std::string msg = "argument with key '" + key + "' not found";
			Error(msg);
		}

		return it->second;
	}

	void AddArgValue(const std::string &key, const std::string &value)
	{
		if (NonKey(value))
		{
			m_args[key].push_back(value);

			// TODO: kostyl'
			const Rule &rule = FindRule(key);

			if (!rule.secondaryKey.empty())
			{
				m_args[rule.secondaryKey].push_back(value);
			}
			//
		}
		else
		{
			std::string msg = "one argument of key '" + key + "' is incorrect";
			Error(msg);
		}
	}

	void CheckDependency(const std::string &key, const Rule &rule)
	{
		std::string dependsOn = rule.dependsOn;

		if (!dependsOn.empty())
		{
			auto it = m_args.find(dependsOn);

			if (it == m_args.end())
			{
				std::string msg = "the argument that argument '"
						+ key + "' depends on is not found";
				Error(msg);
			}
		}
	}

	void Error(const std::string &msg) const
	{
		std::string errMsg = "ArgPP error: ";
		errMsg += msg;
		std::cerr << errMsg << std::endl;
		throw std::exception();
	}

	void AddArgRule(const std::string &key, const Rule &rule)
	{
		if (m_rules.find(key) == m_rules.end())
		{
			m_rules.insert(NamedRule(key, rule));
		}
		else
		{
			std::string msg = "argument with key '" + key + "' already exists";
			Error(msg);
		}
	}

	std::string RetrieveKey(const std::string &rawArg)
	{
		using namespace std;
		string key;
		size_t c = 0;
		bool isOk = false;

		if ((rawArg.size() > 1) && (IsKey(rawArg[c])))
		{
			++c;

			if (IsKey(rawArg[c]))
			{
				isOk = (rawArg.size() > 3) && isalpha(rawArg[++c], locale());
			}
			else
			{
				isOk = (rawArg.size() == 2) && isalpha(rawArg[c], locale());
			}
		}

		if (isOk)
		{
			key = rawArg.substr(c);
		}

		return key;
	}

	void CheckRequiredArgs()
	{
		for (const NamedRule &nrule : m_rules)
		{
			Rule rule = nrule.second;

			if (!rule.isOptional
					&& m_args.find(nrule.first) == m_args.end()
					&& m_args.find(rule.secondaryKey) == m_args.end())
			{
				std::string msg = "required argument with key '"
						+ nrule.first + "' not found";
				Error(msg);
			}
		}
	}

	bool NonKey(const std::string &arg)
	{
		return arg[0] != '-';
	}

	bool IsKey(const char symb)
	{
		return symb == '-';
	}

	template <class T>
	T ConvertTo(const std::string &rawValue) const
	{
		bool ok = true;
		char *end;
		T val;

		if (typeid(T) == typeid(int))
		{
			val = (T)strtol(rawValue.c_str(), &end, 10);
		}
		else if (typeid(T) == typeid(double))
		{
			val = (T)strtod(rawValue.c_str(), &end);
		}
		else
		{
			ok = false;
		}

		if (strlen(end) != 0 || !ok)
		{
			std::string msg = "cannot convert value '" + rawValue
					+ "' to integer";
			Error(msg);
		}

		return val;
	}
};
