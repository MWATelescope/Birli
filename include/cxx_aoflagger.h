#pragma once
#include "rust/cxx.h"
#include <memory>
#include <aoflagger.h>

using namespace std;
using namespace aoflagger;

class CxxImageSet {
friend class CxxAOFlagger;
friend class CxxStrategy;
public:
size_t Width() const;
size_t Height() const;
size_t ImageCount() const;
size_t HorizontalStride() const;
rust::Slice<float> ImageBuffer(size_t imageIndex) const;
private:
CxxImageSet();
CxxImageSet(ImageSet impl);
shared_ptr<ImageSet> pImpl;
};

class CxxFlagMask {
friend class CxxAOFlagger;
friend class CxxStrategy;
public:
size_t Width() const;
size_t Height() const;
size_t HorizontalStride() const;
rust::Slice<bool> Buffer() const;
private:
CxxFlagMask();
CxxFlagMask(FlagMask impl);
shared_ptr<FlagMask> pImpl;
};

class CxxStrategy {
friend class CxxAOFlagger;
public:
unique_ptr<CxxFlagMask> Run(const CxxImageSet& input) const;
unique_ptr<CxxFlagMask> RunExisting(const CxxImageSet& input, const CxxFlagMask& existingFlags) const;
private:
CxxStrategy(Strategy* impl);
mutable Strategy impl;
};

class CxxAOFlagger {
public:
CxxAOFlagger();
void GetVersion(short& major, short& minor, short& subMinor) const;
unique_ptr<CxxImageSet> MakeImageSet(size_t width, size_t height, size_t count, float initialValue, size_t widthCapacity) const;
unique_ptr<CxxFlagMask> MakeFlagMask(size_t width, size_t height, bool initialValue) const;
unique_ptr<CxxStrategy> LoadStrategyFile(const rust::String& filename) const;
rust::String FindStrategyFileGeneric(const rust::String& scenario = "") const;
rust::String FindStrategyFileMWA() const;
private:
// Opaque pointer to aoflagger implementation
shared_ptr<AOFlagger> pImpl;
};

void aoflagger_GetVersion(short& major, short& minor, short& subMinor);
unique_ptr<CxxAOFlagger> cxx_aoflagger_new();


