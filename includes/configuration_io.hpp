#include <string>
#include <sys/stat.h>
#include <memory>
#pragma once

int pushConfig(std::string filename);
int pullConfig(std::string filename);
inline bool fileExists(const std::string &filename)
{
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}
void createIdentifier(std::string &unique_key);
int configureStep();