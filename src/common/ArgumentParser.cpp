#include "ArgumentParser.h"

ArgumentParser::ArgumentParser(ArgPP *parser)
{

}

AngleRange ArgumentParser::GetRange(const Orientation &symmetry)
{
	AngleRange range;

	range.numberZenith = m_parser->GetIntValue("random", 0);
	GetRangeParams(symmetry, "b", range.min.zenith, range.max.zenith);

	range.numberAzimuth = m_parser->GetIntValue("random", 1);
	GetRangeParams(symmetry, "g", range.min.azimuth, range.max.azimuth);

	range.Setup();
	return range;
}

void ArgumentParser::GetRangeParams(const Orientation &symmetry,
									const std::string &key,
									double &min, double &max)
{
	if (m_parser->IsCatched(key))
	{
		min = Orientation::DegToRad(m_parser->GetDoubleValue(key, 0));
		max = Orientation::DegToRad(m_parser->GetDoubleValue(key, 1));
	}
	else
	{
		min = 0;

		if (key == "b")
		{
			max = symmetry.zenith;
		}
		else if (key == "g")
		{
			max = symmetry.azimuth;
		}
		else
		{
			std::cerr << "Error! " << __FUNCTION__;
			throw std::exception();
		}
	}
}

void ArgumentParser::SetArgRules()
{
	int zero = 0;
	m_parser->AddRule("p", '+'); // particle (type, size, ...)
	m_parser->AddRule("ri", 2); // refractive index (Re and Im parts)
	m_parser->AddRule("n", 1); // number of internal reflection
	m_parser->AddRule("pf", zero, true); // particle (filename)
	m_parser->AddRule("fixed", 2, true); // fixed orientarion (beta, gamma)
	m_parser->AddRule("random", 2, true); // random orientarion (beta number, gamma number)
	m_parser->AddRule("go", 0, true); // geometrical optics method
	m_parser->AddRule("po", 0, true); // phisical optics method
	m_parser->AddRule("w", 1, true); // wavelength
	m_parser->AddRule("b", 2, true); // beta range (begin, end)
	m_parser->AddRule("g", 2, true); // gamma range (begin, end)
	m_parser->AddRule("conus", 3, true, "po"); // calculate only backscatter cone (radius, phi, theta)
	m_parser->AddRule("point", zero, true, "po"); // calculate only backscatter poin
	m_parser->AddRule("con20", zero, true, "po");
	m_parser->AddRule("gr", zero, true);
	m_parser->AddRule("tr", 1, true); // file with trajectories
	m_parser->AddRule("all", 0, true); // calculate all trajectories
	m_parser->AddRule("abs", zero, true, "w"); // accounting of absorbtion
	m_parser->AddRule("close", 0, true); // closing of program after calculation
	m_parser->AddRule("o", 1, true); // output folder name
}
