#include "birli/include/cxx_aoflagger.h"
#include "birli/src/lib.rs.h"

using namespace std;
using namespace aoflagger;
using namespace rust;

void aoflagger_GetVersion(short& major, short& minor, short& subMinor)
{
	AOFlagger::GetVersion(major, minor, subMinor);
}

CxxImageSet::CxxImageSet() : pImpl(new ImageSet()) {
}
CxxImageSet::CxxImageSet(ImageSet impl) : pImpl(shared_ptr<ImageSet>(new ImageSet(impl))) {
}
size_t CxxImageSet::Width() const {
	return this->pImpl->Width();
}
size_t CxxImageSet::Height() const {
	return this->pImpl->Height();
}
size_t CxxImageSet::ImageCount() const {
	return this->pImpl->ImageCount();
}
size_t CxxImageSet::HorizontalStride() const {
	return this->pImpl->HorizontalStride();
}
rust::Slice<float> CxxImageSet::ImageBuffer(size_t imageIndex) const {
	rust::Slice<float> slice{this->pImpl->ImageBuffer(imageIndex), Height() * HorizontalStride()};
	return slice;
}

CxxFlagMask::CxxFlagMask() : pImpl(new FlagMask()) {
}
CxxFlagMask::CxxFlagMask(FlagMask impl) : pImpl(shared_ptr<FlagMask>(new FlagMask(impl))) {
}
size_t CxxFlagMask::Width() const {
	return this->pImpl->Width();
}
size_t CxxFlagMask::Height() const {
	return this->pImpl->Height();
}
size_t CxxFlagMask::HorizontalStride() const {
	return this->pImpl->HorizontalStride();
}
rust::Slice<bool> CxxFlagMask::Buffer() const {
	rust::Slice<bool> slice{(bool *)(this->pImpl->Buffer()), Height() * HorizontalStride()};
	return slice;
}

CxxStrategy::CxxStrategy(Strategy* impl) {
    this->impl = std::move(*impl);
}
unique_ptr<CxxFlagMask> CxxStrategy::Run(const CxxImageSet& input) const {
	FlagMask flagmask = this->impl.Run(*(input.pImpl));
	return unique_ptr<CxxFlagMask>(new CxxFlagMask(flagmask));
}
unique_ptr<CxxFlagMask> CxxStrategy::RunExisting(const CxxImageSet& input, const CxxFlagMask& existingFlags) const {
	FlagMask flagmask = this->impl.Run(*(input.pImpl), *(existingFlags.pImpl));
	return unique_ptr<CxxFlagMask>(new CxxFlagMask(flagmask));
}

CxxAOFlagger::CxxAOFlagger() : pImpl(new AOFlagger()) {
}
void CxxAOFlagger::GetVersion(short& major, short& minor, short& subMinor) const {
	this->pImpl->GetVersion(major, minor, subMinor);
}
unique_ptr<CxxImageSet> CxxAOFlagger::MakeImageSet(size_t width, size_t height, size_t count, float initialValue, size_t widthCapacity) const {
	ImageSet imageset = this->pImpl->MakeImageSet(width, height, count, initialValue, widthCapacity);
	return unique_ptr<CxxImageSet>(new CxxImageSet(imageset));
}
unique_ptr<CxxFlagMask> CxxAOFlagger::MakeFlagMask(size_t width, size_t height, bool initialValue) const {
	FlagMask flagmask = this->pImpl->MakeFlagMask(width, height, initialValue);
	return unique_ptr<CxxFlagMask>(new CxxFlagMask(flagmask));
}
rust::String CxxAOFlagger::FindStrategyFileGeneric(const rust::String& scenario) const {
	return this->pImpl->FindStrategyFile(TelescopeId::GENERIC_TELESCOPE, std::string(scenario));
}
rust::String CxxAOFlagger::FindStrategyFileMWA() const {
	return this->pImpl->FindStrategyFile(TelescopeId::MWA_TELESCOPE);
}
unique_ptr<CxxStrategy> CxxAOFlagger::LoadStrategyFile(const rust::String& filename) const {
	Strategy strategy = this->pImpl->LoadStrategyFile(std::string(filename));
	return unique_ptr<CxxStrategy>(new CxxStrategy(&strategy));
}

unique_ptr<CxxAOFlagger> cxx_aoflagger_new() {
	return unique_ptr<CxxAOFlagger>(new CxxAOFlagger());
};

