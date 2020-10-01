#include "FileHandler.h"

#include <iostream>
#include <fstream>


std::string FileHandler::getFileExtension(const std::string& sstrFileName)
{
	if(sstrFileName.empty())
	{
		return "";
	}

	const size_t pos = sstrFileName.rfind(".");
	if(pos != std::string::npos)
	{
		return sstrFileName.substr(pos+1);
	}

	return "";
}

bool FileHandler::fileExist(const std::string& sstrFileName)
{
	if(sstrFileName.empty())
	{
		return false;
	}

	std::fstream inStream;
	inStream.open(sstrFileName, std::ios::in);

	if(inStream.is_open())
	{
		inStream.close();
		return true;
	}
	
	return false;
}

std::string FileHandler::getFileName(const std::string& sstrFullFileName)
{
	std::string sstrFileName = "";

	if(sstrFullFileName.empty())
	{
		return sstrFileName;
	}

	if(!FileHandler::fileExist(sstrFullFileName))
	{
		return sstrFileName;
	}

	size_t posEnd = sstrFullFileName.rfind(".");
	if(posEnd == std::string::npos)
	{
		posEnd = sstrFullFileName.size();
	}

	size_t posStart = sstrFullFileName.rfind("/");
	if(posStart != std::string::npos)
	{
		sstrFileName = sstrFullFileName.substr(posStart+1, posEnd-posStart-1);
		return sstrFileName;
	}

	posStart = sstrFullFileName.rfind("\\");
	if(posStart != std::string::npos)
	{
		sstrFileName = sstrFullFileName.substr(posStart+1, posEnd-posStart-1);
		return sstrFileName;
	}

	return sstrFileName;
}

std::string FileHandler::getFilePath(const std::string& sstrFileName)
{
	std::string sstrFilePath = "";

	if(sstrFileName.empty())
	{
		return sstrFilePath;
	}

	if(!FileHandler::fileExist(sstrFileName))
	{
		return sstrFilePath;
	}

	size_t pos = sstrFileName.rfind("/");
	if(pos != std::string::npos)
	{
		sstrFilePath = sstrFileName.substr(0, pos);
		return sstrFilePath;
	}

	pos = sstrFileName.rfind("\\");
	if(pos != std::string::npos)
	{
		sstrFilePath = sstrFileName.substr(0, pos);
		return sstrFilePath;
	}

	return sstrFilePath;
}


