/** @file aoflagger.h @brief Main AOFlagger header file.
 * @author André Offringa offringa@gmail.com
 * @copyright by A.R. Offringa under the GPL version 3
 */

#ifndef AOFLAGGER_INTERFACE_H
#define AOFLAGGER_INTERFACE_H

#include <cstring>
#include <string>
#include <memory>

/** @brief Contains all the public types used by the AOFlagger.
 *
 * See the @ref AOFlagger class description for details.
 * @author André Offringa offringa@gmail.com
 */
namespace aoflagger {

/** @brief Strategy identifier for the supported telescopes.
 *
 * If you have an optimized strategy for an unlisted telescope, please
 * contact me.
 * @sa AOFlagger::FindStrategyFile().
 * @since Version 3.0
 */
enum class TelescopeId {
	/** @brief Most generic strategy. */
	GENERIC_TELESCOPE,
	/** @brief The AARTFAAC telescope, correlating the superterp antennas of LOFAR. */
	AARTFAAC_TELESCOPE,
	/** @brief The WSRT telescope with the Apertif focal plane array receiver system. */
	APERTIF_TELESCOPE,
	/** @brief Arecibo radio telescope, the 305 m telescope in Puerto Rico. */
	ARECIBO_TELESCOPE,
	/** @brief Bighorns, instrument aimed at achieving an averaged all-sky measurement of the Epoch of Reionisation signal. */
	BIGHORNS_TELESCOPE,
	/** @brief JVLA, the Jansky Very Large Array in New Mexico. */
	JVLA_TELESCOPE,
	/** @brief LOFAR. the Low-Frequency Array in Europe. */
	LOFAR_TELESCOPE,
	/** @brief MWA, the Murchison Widefield Array in Western Australia. */
	MWA_TELESCOPE,
	/** @brief Parkes, the single dish telescope in New South Wales. */
	PARKES_TELESCOPE,
	/** @brief WSRT, the Westerbork Synthesis Radio Telescope in the Netherlands. */
	WSRT_TELESCOPE
};

/** @brief A set of time-frequency 'images' which together contain data for one
 * correlated baseline or dish.
 *
 * The class either holds 1, 2, 4 or 8 images. These images have time on the
 * x-axis (most rapidly changing index) and frequency on the y-axis. The
 * cells specify flux levels, which do not need to have been calibrated.
 *
 * If the set contains only one image, it specifies amplitudes of a single
 * polarization. If it contains two images, it specifies the real and imaginary
 * parts of a single polarization. With four images, it contains the real
 * and imaginary values of two polarizations (ordered real pol A, imag pol A,
 * real pol B, imag pol B). With eight images, it contains complex values for
 * four correlated polarizations (ordered real pol A, imag pol A, real pol B,
 * ... etc).
 *
 * @note When accesses the image data, note that there might be more items on one row
 * than the width of the image. The rows are padded to align them e.g. for
 * SSE instructions. Use @ref HorizontalStride() to get the actual number of
 * floats per row.
 */
class ImageSet
{
public:
	friend class AOFlagger;
	friend class QualityStatistics;
	friend class Strategy;

	/** @brief Construct an empty ImageSet.
	 *
	 * The only operations allowed on an empty ImageSet are to assign to it. Use
	 * AOFlagger::MakeImageSet() to construct a non-empty ImageSet.
	 */
	ImageSet();

	/** @brief Copy the image set. Only references to images are copied. */
	ImageSet(const ImageSet& sourceImageSet);

	/** @brief Move from the image set.
	 * @since Version 2.10
	 */
	ImageSet(ImageSet&& sourceImageSet);

	/** @brief Destruct image set. Destroys its images if no longer referenced. */
	~ImageSet();

	/** @brief Assign to this image set. Only references to images are copied. */
	ImageSet& operator=(const ImageSet& sourceImageSet);

	/** @brief Move assign to this image set.
	 * @since Version 2.10
	 */
	ImageSet& operator=(ImageSet&& sourceImageSet);

	/** @brief Get access to the data buffer of an image.
	 * @param imageIndex Index of image. See class description for ordering.
	 * \note Rows are padded, see @ref HorizontalStride().
	 */
	float* ImageBuffer(size_t imageIndex);

	/** @brief Get constant access to the data buffer of an image.
	 * @param imageIndex Index of image. See class description for ordering.
	 * \note Rows are padded, see @ref HorizontalStride().
	 */
	const float* ImageBuffer(size_t imageIndex) const;

	/** @brief Get width (number of time steps) of images. */
	size_t Width() const;

	/** @brief Get height (number of frequency channels) of images. */
	size_t Height() const;

	/** @brief Get number of images, see class description for details. */
	size_t ImageCount() const;

	/** @brief Get total number of floats in one row.
	 *
	 * Row might have been padded to allow for
	 * SSE instructions and other optimizations. Therefore, one should
	 * add the horizontal stride to a data pointer to get the float in the next
	 * row (channel).
	 *
	 * Example:
	 * @code{.cpp}
	 *(ImageSet::ImageBuffer(imageIndex) + x + y * ImageSet::HorizontalStride())
	 * @endcode
	 * will return the value at position x,y.
	 */
	size_t HorizontalStride() const;

	/** @brief Set all samples to the specified value.
	 * @param newValue The new value for all values of all images in the set.
	 * @since 2.5
	 */
	void Set(float newValue);

	/** @brief Resize the image without reallocating new memory.
	 *
	 * This function allows to quickly change the dimension of the images in the
	 * imageset. The new width has to fit in the image capacity as specified
	 * during creation. When flagging many images of "almost" the same size, using
	 * this method to change the size of images is drastically faster compared
	 * to freeing and then allocating new images. It was added after rather
	 * severe memory fragmentation problems in the Cotter MWA pipeline.
	 * @param newWidth The new width of the images. Should satisfy newWidth <= HorizontalStride().
	 * @since 2.5
	 */
	void ResizeWithoutReallocation(size_t newWidth) const;

private:
	ImageSet(size_t width, size_t height, size_t count);

	ImageSet(size_t width, size_t height, size_t count, float initialValue);

	ImageSet(size_t width, size_t height, size_t count, size_t widthCapacity);

	ImageSet(size_t width, size_t height, size_t count, float initialValue, size_t widthCapacity);

	static void assertValidCount(size_t count);

	std::unique_ptr<class ImageSetData> _data;
};

/** @brief A two-dimensional flag mask.
 *
 * The flag mask specifies which values in an @ref ImageSet are flagged.
 * A value @c true means a value is flagged, i.e., contains RFI and should
 * not be used in further data processing (calibration, imaging, etc.).
 * A flag denotes that the value at that time-frequency position should
 * be ignored for all polarizations. This normally makes sense, because if one
 * polarization is contaminated by RFI, all polarizations are probably
 * affected. Also, solving for Stokes matrices during calibration might
 * not work well when the polarizations are not flagged equally.
 *
 * If polarization-specific flags are needed, one could run the flagger on
 * each polarization individually. However, note that some algorithms, like
 * the morphological scale-invariant rank operator (SIR operator), work best
 * when seeing the flags from all polarizations.
 *
 * @note When accesses the flag data, note that there might be more items on one row
 * than the width of the mask. The rows are padded to align them e.g. for
 * SSE instructions. Use @ref HorizontalStride() to get the actual number of
 * bools per row.
 */
class FlagMask
{
public:
	friend class AOFlagger;
	friend class QualityStatistics;
	friend class Strategy;

	/** @brief Construct an empty FlagMask.
	 * The properties of an empty FlagMask can not be accessed.
	 */
	FlagMask();

	/** @brief Copy a flag mask. Only copies a reference, not the data. */
	FlagMask(const FlagMask& sourceMask);

	/** @brief Move construct a flag mask.
	 * @since Version 2.10
	 */
	FlagMask(FlagMask&& sourceMask);

	/** @brief Copy assignment.
	 * @since Version 2.10
	 */
	FlagMask& operator=(const FlagMask& source);

	/** @brief Move assignment.
	 * @since Version 2.10
	 */
	FlagMask& operator=(FlagMask&& source);

	/** @brief Destroy a flag mask. Destroys mask data if no longer references. */
	~FlagMask();

	/** @brief Get the width of the mask. */
	size_t Width() const;

	/** @brief Get the height of the mask. */
	size_t Height() const;

	/** @brief Get total number of bools in one row.
	 *
	 * Row might have been padded to allow for
	 * SSE instructions and other optimizations. Therefore, one should
	 * add the horizontal stride to a data pointer to get the flags in
	 * the next row (channel).
	 *
	 * Example:
	 * @code{.cpp}
	 *(FlagMask::Buffer() + x + y * Buffer::HorizontalStride())
	 * @endcode
	 * will return the flag value at position x,y.
	 */
	size_t HorizontalStride() const;

	/** @brief Get access to the data buffer.
	 * @note The buffer is padded, see @ref HorizontalStride(). */
	bool* Buffer();

	/** @brief Get constant access to the data buffer.
	 * @note The buffer is padded, see @ref HorizontalStride(). */
	const bool* Buffer() const;

private:
	FlagMask(size_t width, size_t height);
	FlagMask(size_t width, size_t height, bool initialValue);

	std::unique_ptr<class FlagMaskData> _data;
};

/** @brief Holds a flagging strategy.
 *
 * Default "stock" strategies can be found with
 * @ref AOFlagger::FindStrategyFile(), and these or custom Lua files can
 * be loaded from disc with @ref AOFlagger::LoadStrategyFile().
 * A user can create strategies with the @c rfigui tool that is part
 * of the aoflagger package.
 *
 * When flagging a large number of baselines it is recommended to use multiple threads.
 * This class is itself not thread save, but it is safe to use different Strategy objects
 * from different thread contexts.
 */
class Strategy
{
public:
	friend class AOFlagger;

	/** @brief Construct an empty strategy.
	 *
	 * The only operations allowed on an empty Strategy are to assign to it. Use
	 * e.g. AOFlagger::LoadStrategyFile() to construct a non-empty Strategy.
	 * @since Version 3
	 */
	Strategy();

	/** @brief Move construct a strategy.
	 * @since Version 2.10
	 */
	Strategy(Strategy&& sourceStrategy);

	/** @brief Destruct strategy. */
	~Strategy();

	/** @brief Move assign to strategy.
	 * @since Version 2.10
	 */
	Strategy& operator=(Strategy&& sourceStrategy);

	/** @brief Run the flagging strategy on the given data.
	 *
	 * The Lua strategy is executed single-threaded. This function is not
	 * thread safe: To flag multiple imagesets simultaneously from different threads,
	 * it is necessary to create a Strategy object for each thread.
	 *
	 * @param input The data to run the flagger on.
	 * @return The flags identifying bad (RFI contaminated) data.
	 * @since 3.0
	 */
	FlagMask Run(const ImageSet& input);

	/** @brief Run the flagging strategy on the given data with existing flags.
	 *
	 * This method is similar to @ref Run(const ImageSet&), except
	 * that it will pass existing flags (e.g. as set by the correlator)
	 * to the flagging strategy, which in the case of bad data can do a better
	 * job of finding RFI in the good data.
	 * @p input parameter. The @p strategy parameter can be the
	 * same for different threads.
	 * @param input The data to run the flagger on.
	 * @param existingFlags Flags that indicate what data are bad.
	 * @return A flag mask that identifies bad (RFI contaminated) data.
	 * @since 3.0
	 */
	FlagMask Run(const ImageSet& input, const FlagMask& existingFlags);

private:
	Strategy(const std::string& filename, class AOFlagger* aoflagger);
	Strategy(const Strategy& sourceStrategy) = delete;
	Strategy& operator=(const Strategy& sourceStrategy) = delete;

	FlagMask run(const ImageSet& input, const FlagMask* existingFlags);

	std::unique_ptr<class StrategyData> _data;
	class AOFlagger* _aoflagger;
};

/** @brief Statistics that can be collected online and saved to a measurement set.
 *
 * It is useful to collect some statistics during flagging, because all data goes through
 * memory at highest resolution. This class contains the collected statistics and
 * some meta data required for collecting. It can be created with
 * @ref AOFlagger::MakeQualityStatistics(). Statistics can be added to it with
 * @ref CollectStatistics(), and saved to disk with
 * @ref WriteStatistics().
 *
 * This class does not allow viewing or modifying statistics, it only contains the most
 * basic form to collect statistics during flagging and writing them in the (well-defined)
 * quality statistic tables format. These statistics can be viewed interactively with
 * the @c aoqplot tool.
 *
 * Collecting statistics is not as expensive as flagging, but still takes some time, so it
 * is recommended to use multiple threads for collecting as well. This class is however not
 * thread save, but it is okay to use different QualityStatistics objects from different
 * thread contexts. During finalization, the different objects can be combined with the
 * operator+=() method, and then in full written to the measurement set.
 */
class QualityStatistics
{
public:
	friend class AOFlagger;

	/** Construct a QualityStatistics with null state.
	 *
	 * An object created by this constructor can only be assigned to.
	 * @since Version 2.13
	 */
	QualityStatistics();

	/** @brief Copy the object. This is fast; only references are copied. */
	QualityStatistics(const QualityStatistics& sourceQS);

	/** @brief Move construct the object.
	 * @since Version 2.10
	 */
	QualityStatistics(QualityStatistics&& sourceQS);

	/** @brief Destruct the object. Data is destroyed if no more references exist. */
	~QualityStatistics();

	/** @brief Assign to this object. This is fast; only references are copied. */
	QualityStatistics& operator=(const QualityStatistics& sourceQS);

	/** @brief Move-assign this object. This is fast; only references are moved.
	 * @since Version 2.10
	 */
	QualityStatistics& operator=(QualityStatistics&& sourceQS);

	/** @brief Combine the statistics from the given object with the statistics in this object.
	 *
	 * This is a relative expensive operation, so should only be used scarsely. It can be used
	 * to combine the results of different threads, as explained in the class description.
	 *
	 * It is safe to combine quality statistics with different meta data (scan time count, channel
	 * count, etc.). When using this object again during collecting (see @ref CollectStatistics()),
	 * after combining it with another object, it will still use the meta data it was initialized with.
	 */
	QualityStatistics& operator+=(const QualityStatistics& rhs);

	/** @brief Collect statistics from time-frequency images and masks.
	 *
	 * This will update the statistics in this object so that it
	 * represents the combination of previous collected data and the newly
	 * given data.
	 *
	 * This function can be called from different thread context, as long as each
	 * thread uses its own QualityStatistics object.
	 * See the @ref QualityStatistics class documentation for further multithreading info.
	 * @param imageSet Data to collect statistics from
	 * @param rfiFlags Flags set by the automatic RFI detector
	 * @param correlatorFlags Flags that were set prior to RFI detector, e.g. because of
	 * a broken antenna or correlator hickup.
	 * @param antenna1 Index of the first antenna involved in this baseline.
	 * @param antenna2 Index of the second antenna involved in this baseline.
	 * @since 3.0
	 */
	void CollectStatistics(const ImageSet& imageSet, const FlagMask& rfiFlags, const FlagMask& correlatorFlags, size_t antenna1, size_t antenna2);

	/** @brief Write collected statistics in standard tables to a measurement set.
	 * @param measurementSetPath Path to measurement set to which the statistics will
	 * be written.
	 * @since 3.0
	 */
	void WriteStatistics(const std::string& measurementSetPath) const;

private:
	QualityStatistics(const double* scanTimes, size_t nScans, const double* channelFrequencies, size_t nChannels, size_t nPolarizations, bool computeHistograms);

	std::unique_ptr<class QualityStatisticsData> _data;
};

/**
 * @brief A base class which callers can inherit from to be able to receive
 * progress updates and error messages.
 *
 * A status listener should be thread safe when the Run() method is called in parallel
 * with the same StatusListener object.
 */
class StatusListener
{
public:
	/**
	 * @brief Virtual destructor.
	 */
	virtual ~StatusListener() { }
	/**
	 * @brief This virtual method is called when a new task is started.
	 *
	 * Typically, a client could display a message saying that the given task 'description' is started.
	 * @param description Description of the task, e.g. "SumThreshold".
	 * @since 3.0
	 */
	virtual void OnStartTask(const std::string& description)
	{ }
	/**
	 * @brief Called to update current progress.
	 *
	 * This can be used to display a progress bar if the strategy would take a lot of time.
	 * @param progress Current progress
	 * @param maxProgress Progress that is required to finish the current task.
	 */
	virtual void OnProgress(size_t progress, size_t maxProgress)
	{ }
	/**
	 * @brief Called when detection has completely finished.
	 *
	 * This function will always be called exactly once on success. It is not called when an
	 * exception occurred. @ref OnException() is called in that case.
	 */
	virtual void OnFinish()
	{ }
	/**
	 * @brief Called when an exception occurs during execution of the strategy.
	 *
	 * This can occur when the Lua script throws an error.
	 * @param thrownException The exception that was thrown.
	 */
	virtual void OnException(std::exception& thrownException) = 0;
};

/** @brief Main class for access to the flagger functionality.
 *
 * Software using the flagger should first create an instance of the @ref AOFlagger
 * class, from which other actions can be initiated.
 *
 * ### Overview
 *
 * To flag a data set:
 * - Create the AOFlagger instance
 * - To use a stock strategy, call FindStrategyFile()
 * - Load and parse the strategy with LoadStrategyFile(), once for each thread that will
 * run.
 * - Create data buffers with MakeImageSet()
 * - For each correlated baseline or dish:
 * - - Fill the images with data from this correlated baseline or dish
 * - - Call Strategy::Run() with the created Strategy and ImageSet
 * - - Process the data that was returned in the FlagMask.
 *
 * Optionally, it is possible to assemble quality statistics that can be written to
 * the measurement set in the standard format that e.g. the @c aoqplot tool can read.
 * To do this:
 * - Create (once) a quality statistics object with MakeQualityStatistics().
 * - After flagging a baseline, add it to the statistics object with
 *   QualityStatistics::CollectStatistics().
 * A "correlator mask" can be specified that describes which flags are not due
 * to RFI but caused by different things.
 * - When a full set is processed, store the statistics with WriteStatistics().
 *
 * To flag multiple baselines, the Strategy and/or ImageSet objects can be reused.
 *
 * ### Thread safety
 *
 * Each Strategy object runs in its own context, and will not perform any unsynchronised
 * writes to global variables. It is therefore safe to call Strategy::Run() from
 * different threads, as long as each thread uses its own Strategy and ImageSet instances.
 * QualityStatistics::CollectStatistics() is also thread safe, as long as different QualityStatistics
 * instances are passed. For multi-threading, each thread should collect into
 * its own QualityStatistics object. When finished, these can be combined with
 * QualityStatistics::operator+=().
 *
 * It is safe to create multiple AOFlagger instances, but not recommended.
 *
 * ### Data order
 *
 * A common problem for integrating the flagger, is that data are stored in a
 * different order: the time dimension
 * is often the direction with the slowest increasing indices. Because the flagger
 * needs one baseline at a time, this requires reordering the data. As long as the
 * data fits in memory, this reordering is quite straightforward. When this is not the
 * case, the data could be split into sub-bands and/or time windows.
 * Next, these parts can be passed to the flagger and recombined later (if desired).
 *
 * To decide how to split, keep in mind that the flagger
 * works best when both a lot of channels and a lot of
 * timesteps are available. As an example: LOFAR splits into subbands of ~64 channels, and
 * the default processing with NDPPP loads as many timesteps as possible
 * in memory for flagging. Typically, this means at least a
 * few hundred of timesteps are processed at a time (with 1-3s per timestep), and
 * this seems to work well.
 *
 * The 'aoflagger' executable flags by default on the full measurement set.
 * For sets that are larger than memory, a mode is used in
 * which the data is reordered to disk before the actual flagging starts. It turns
 * out that this is much faster than reading each baseline directly from the set, and
 * it simultaneously produces the best flagging accuracy, so
 * if enough processing power is available to do so, this is another approach.
 */
class AOFlagger
{
public:
	/** @brief Create and initialize the flagger main class. */
	AOFlagger() : _statusListener(nullptr) { }

	/** @brief Destructor. */
	~AOFlagger() { }

	/** @brief Create a new uninitialized @ref ImageSet with specified specs.
	 *
	 * The float values will not be initialized.
	 * @param width Number of time steps in images
	 * @param height Number of frequency channels in images
	 * @param count Number of images in set (see class description
	 * of @ref ImageSet for image order).
	 * @return A new ImageSet.
	 */
	ImageSet MakeImageSet(size_t width, size_t height, size_t count)
	{
		return ImageSet(width, height, count);
	}

	/** @brief Create a new uninitialized @ref ImageSet with specified specs.
	 *
	 * The float values will not be initialized.
	 * @param width Number of time steps in images
	 * @param height Number of frequency channels in images
	 * @param count Number of images in set (see class description
	 * of @ref ImageSet for image order).
	 * @param widthCapacity Allow for enlarging image to this size, @sa ImageSet::ResizeWithoutReallocation()
	 * @return A new ImageSet.
	 * @since 2.6
	 */
	ImageSet MakeImageSet(size_t width, size_t height, size_t count, size_t widthCapacity)
	{
		return ImageSet(width, height, count, widthCapacity);
	}

	/** @brief Create a new initialized @ref ImageSet with specified specs.
	 * @param width Number of time steps in images
	 * @param height Number of frequency channels in images
	 * @param count Number of images in set (see class description
	 * of @ref ImageSet for image order).
	 * @param initialValue Initialize all pixels with this value.
	 * @return A new ImageSet.
	 */
	ImageSet MakeImageSet(size_t width, size_t height, size_t count, float initialValue)
	{
		return ImageSet(width, height, count, initialValue);
	}

	/** @brief Create a new initialized @ref ImageSet with specified specs.
	 * @param width Number of time steps in images
	 * @param height Number of frequency channels in images
	 * @param count Number of images in set (see class description
	 * of @ref ImageSet for image order).
	 * @param initialValue Initialize all pixels with this value.
	 * @param widthCapacity Allow for enlarging image to this size, @sa ImageSet::ResizeWithoutReallocation()
	 * @return A new ImageSet.
	 * @since 2.6
	 */
	ImageSet MakeImageSet(size_t width, size_t height, size_t count, float initialValue, size_t widthCapacity)
	{
		return ImageSet(width, height, count, initialValue, widthCapacity);
	}

	/** @brief Create a new uninitialized @ref FlagMask with specified dimensions.
	 * @param width Width of mask (number of timesteps)
	 * @param height Height of mask (number of frequency channels)
	 * @return A new FlagMask.
	 */
	FlagMask MakeFlagMask(size_t width, size_t height)
	{
		return FlagMask(width, height);
	}

	/** @brief Create a new initialized @ref FlagMask with specified dimensions.
	 * @param width Width of mask (number of timesteps)
	 * @param height Height of mask (number of frequency channels)
	 * @param initialValue Value to initialize the mask to.
	 * @return A new FlagMask.
	 */
	FlagMask MakeFlagMask(size_t width, size_t height, bool initialValue)
	{
		return FlagMask(width, height, initialValue);
	}

	/** @brief Find a Lua strategy for a specific telescope.
	 *
	 * The scenario name can be used to
	 * distinguish different strategies for the same telescope, for example
	 * to distinguish between different bands, different processing stages (e.g.
	 * before or after averaging), different antennas (LOFAR HBA vs LBA), etc.
	 *
	 * This will search the data directory of the installation path of
	 * aoflagger for a file named &lt;Telescope_Name&gt;-&lt;scenario&gt;.lua, or
	 * &lt;Telescope_Name&gt;-default.lua in case scenario is empty.
	 *
	 * @param telescopeId Identifies the telescope to optimize the strategy for.
	 * @param scenario A scenario name that should be searched for.
	 * @returns Filename of strategy, or empty string if none found.
	 * @since Version 3.0
	 */
	std::string FindStrategyFile(enum TelescopeId telescopeId=TelescopeId::GENERIC_TELESCOPE, const std::string& scenario = "");

	/** @brief Load a strategy from disk.
	 *
	 * The best way to create strategies is to use the @c rfigui tool. In case you have optimized
	 * strategies for an unlisted telescope or for new scenarios, please consider
	 * providing the file so I can add them to the repository.
	 *
	 * @param filename Full pathname to .lua strategy file.
	 * @return The new @ref Strategy.
	 * @since Version 3.0
	 */
	Strategy LoadStrategyFile(const std::string& filename)
	{
		return Strategy(filename, this);
	}

	/**
	 * @brief Create a new object for collecting statistics.
	 * @param scanTimes Array with times. The number of elements should match the
	 * dimension of the time axis in calls to QualityStatistics::CollectStatistics().
	 * Each time is a MJD time in seconds.
	 * @param nScans Number of elements in the @c scanTimes array.
	 * @param channelFrequencies Frequency in Hz of each channel. The number of elements
	 * should match the frequency axis in calls to QualityStatistics::CollectStatistics().
	 * @param nChannels Number of elements in the @c channelFrequencies array.
	 * @param nPolarizations Number of polarizations in the set (1, 2 or 4).
	 *
	 * See the QualityStatistics class description for info on multithreading and/or combining statistics
	 * with different meta data. The meta data that is passed to this method will be used for all
	 * calls to QualityStatistics::CollectStatistics(). No histograms will be computed.
	 */
	QualityStatistics MakeQualityStatistics(const double* scanTimes, size_t nScans, const double* channelFrequencies, size_t nChannels, size_t nPolarizations);

	/** @brief Create a new object for collecting statistics, possibly with histograms.
	 *
	 * See other overload of MakeQualityStatistics() for info.
	 * @since Version 2.6
	 */
	QualityStatistics MakeQualityStatistics(const double* scanTimes, size_t nScans, const double* channelFrequencies, size_t nChannels, size_t nPolarizations, bool computeHistograms);

	/** @brief Get the AOFlagger version number as a string.
	 * @returns The version number, formatted like '1.2.3-subtitle', e.g. '3.0-alpha'.
	 * @since Version 2.6
	 */
	static std::string GetVersionString();

	/** @brief Get the AOFlagger version number separated in major, minor and subminor fields.
	 * @param major Most significant number of the version, e.g. '1' for version '1.2.3'. This
	 * number is only incremented in major changes of the flagger.
	 * @param minor Minor number of the version, e.g. '2' for version '1.2.3'. This number
	 * is incremented for every public release.
	 * @param subMinor Subminor number of the version, e.g. '3' for version '1.2.3', or zero if
	 * the current version has no subminor number. This number is incremented for internal releases
	 * or small bug fixes.
	 * @since Version 2.6
	 */
	static void GetVersion(short& major, short& minor, short& subMinor);

	/** @brief Get the date this version was released as a string.
	 * @returns The version date formatted like "1982-05-08".
	 * @since Version 2.6
	 */
	static std::string GetVersionDate();

	/**
	 * @brief Set a handler for progress updates and exceptions.
	 *
	 * By default, exceptions will be reported to stderr and progress updates
	 * will be ignored. If an application needs to handle either of these
	 * themselves, they can override a @ref StatusListener that handles these
	 * events and call this method to enable receiving the events.
	 * This method is not thread safe.
	 * @param statusListener The handler that will receive the status updates.
	 * @since Version 2.7
	 */
	void SetStatusListener(StatusListener* statusListener)
	{
		_statusListener = statusListener;
	}

private:
	friend class Strategy;

	/** @brief It is not allowed to copy this class
	 */
	AOFlagger(const AOFlagger&) = delete;

	/** @brief It is not allowed to assign to this class
	 */
	void operator=(const AOFlagger&) = delete;

	StatusListener* _statusListener;
};

} // namespace aoflagger

#endif
