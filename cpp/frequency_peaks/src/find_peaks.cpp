#include <cstdint>
#include <cmath>
#include <cassert>

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>

#if defined(_MSC_VER)
#	define DLL_API extern "C" __declspec(dllexport)
#	define CDECL __cdecl
#elif defined(__GNUC__)
#	define DLL_API extern "C" __attribute__((visibility("default")))
#	define CDECL __attribute__((__cdecl__))
#else
#	define DLL_API extern "C"
#	define CDECL
#	warning Unhandled dll export attributes
#endif

using iter = std::vector<double>::const_iterator;

static iter find_first_peak_impl(iter it, iter const end, size_t avg_len)
{
	assert(end >= it);
	assert((size_t)(end - it) >= 3 * avg_len);

	iter trail = it;
	iter middle = it;
	double trail_sum = 0.0;
	for (size_t i = 0; i < avg_len; ++i, ++middle)
	{
		assert(middle != end);
		trail_sum += *middle;
	}

	iter advance = middle;
	double middle_sum = 0.0;
	for (size_t i = 0; i < avg_len; ++i, ++advance)
	{
		assert(advance != end);
		middle_sum += *advance;
	}

	iter front = advance;
	double advance_sum = 0.0;
	for (size_t i = 0; i < avg_len; ++i, ++front)
	{
		assert(front != end);
		advance_sum += *front;
	}

#define next()               \
do                           \
{                            \
    trail_sum -= *trail;     \
    ++trail;                 \
    trail_sum += *middle;    \
    middle_sum -= *middle;   \
    ++middle;                \
    middle_sum += *advance;  \
    advance_sum -= *advance; \
    ++advance;               \
    advance_sum += *front;   \
    ++front;                 \
} while(false)

	while (true)
	{
		double threshold = 6.0 * avg_len;
		if (
			middle_sum > trail_sum + threshold
			&& middle_sum > advance_sum + threshold
		)
		{
			return std::max_element(trail, front);
		}

		if (front == end)
		{
			break;
		}
		next();
	}

	return end;
}

static iter find_first_peak(iter begin, iter const end, double df)
{
	std::vector<iter> peaks;

	auto const range_begin = std::max(static_cast<size_t>(2  / df), 1ull);
	auto const range_end   = std::max(static_cast<size_t>(20 / df), 2ull);

	begin += static_cast<int64_t>(52.0 / df);

	for (size_t i = range_begin; i < range_end; ++i)
	{
		peaks.push_back(find_first_peak_impl(begin, end, i));
	}

//	std::cout << std::boolalpha;
//	for (size_t i = 0; i < peaks.size(); ++i)
//	{
//		std::cout << "peaks[" << i << "] (" << i + 2 << ") = " << peaks[i]._Ptr << "; end? " << (peaks[i] == end);
//		if (peaks[i] != end)
//		{
//			std::cout << ' ' << (peaks[i] - begin) << "Hz";
//		}
//		std::cout << '\n';
//	}

	size_t const count = end - begin;
	std::vector<size_t> possible_peaks(count, 0);

	for (size_t i = 0; i < peaks.size(); ++i)
	{
		if (peaks[i] != end)
		{
			++possible_peaks[peaks[i] - begin];
		}
	}

	auto const most_common = std::max_element(possible_peaks.begin(), possible_peaks.end());
	return (most_common - possible_peaks.begin()) + begin;
}

static std::vector<size_t> find_peaks(std::vector<double> const &values, size_t max_peaks, double df)
{
	std::vector<size_t> result;

	auto const values_begin = values.begin();
	auto const values_end   = values.end();

	auto const first_peak = find_first_peak(values_begin, values_end, df);
	if (first_peak != values_end)
	{
		result.push_back(first_peak - values_begin);
	}

	auto prev_peak_it    = first_peak;
	auto const first_peak_index = first_peak - values_begin;

	for (size_t i = 1; i < max_peaks; ++i)
	{
		auto const begin = prev_peak_it + (first_peak_index / 2);
		auto const end   = begin + first_peak_index;

		auto it = std::max_element(begin, end);
		result.push_back(it - values_begin);
		prev_peak_it = it;
	}

	return result;
}

DLL_API void CDECL get_peaks(
	double in[], size_t in_size,
	double out_freq[], size_t out_freq_max_size,
	double out_amp[], size_t out_amp_max_size,
	double df
)
{
	auto const out_size = std::min(out_amp_max_size, out_freq_max_size);

	if (out_size == 0)
	{
		return;
	}

	std::vector<double> values;
	values.resize(in_size);

	std::transform(in, in + in_size, values.begin(), [](double val) {
		return 20 * std::log10(val);
	});

	auto const peaks = find_peaks(values, out_size, df);

	for (size_t i = 0; i < out_size && i < peaks.size(); ++i)
	{
		out_freq[i] = peaks[i] * df;
		out_amp[i]  = in[peaks[i]];
	}
}
