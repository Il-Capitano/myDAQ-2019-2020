#ifdef __GNUC__
#pragma GCC visibility push(hidden)
#endif // gcc

#include <cstdint>
#include <cmath>
#include <cassert>

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <fstream>

#if defined(_MSC_VER) // MSVC
#	define DLL_API extern "C" __declspec(dllexport)
#	define CDECL __cdecl
#elif defined(__GNUC__) // gcc
#	define DLL_API extern "C" __attribute__((visibility("default")))
#	define CDECL __attribute__((__cdecl__))
#else
#	define DLL_API extern "C"
#	define CDECL
#	warning Unhandled dll export attributes
#endif

using iter = std::vector<double>::const_iterator;

static iter find_first_peak_impl(iter it, iter const end, size_t avg_len, double search_threshold)
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
		double threshold = search_threshold * avg_len;
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

static iter find_first_peak(iter begin, iter const end, double df, double threshold, double min_search_freq, double max_search_freq)
{
	std::vector<iter> peaks;

	auto const range_begin = std::max(static_cast<size_t>(min_search_freq / df), 1ull);
	auto const range_end   = std::max(static_cast<size_t>(max_search_freq / df), 2ull);

	begin += static_cast<int64_t>(52.0 / df);

	for (size_t i = range_begin; i < range_end; ++i)
	{
		peaks.push_back(find_first_peak_impl(begin, end, i, threshold));
	}

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

static size_t get_half_value_width(iter begin, iter end, iter max_it)
{
	iter left = begin;
	iter right = end;

	for (auto it = max_it; it != begin;)
	{
		--it;
		if (*it <= *max_it - 20.0 * std::log10(2))
		{
			left = it;
			break;
		}
	}

	for (auto it = max_it + 1; it != end;)
	{
		if (*it > *max_it - 20.0 * std::log10(2))
		{
			right = it;
			break;
		}
		++it;
	}

	return right - left;
}

static std::pair<std::vector<size_t>, std::vector<size_t>> find_peaks(std::vector<double> const &values, size_t max_peaks, double df, double threshold, double min_search_freq, double max_search_freq)
{
	std::pair<std::vector<size_t>, std::vector<size_t>> result;
	auto &peaks  = result.first;
	auto &errors = result.second;

	auto const values_begin = values.begin();
	auto const values_end   = values.end();

	auto const first_peak = find_first_peak(values_begin, values_end, df, threshold, min_search_freq, max_search_freq);
	auto const first_peak_index = first_peak - values_begin;
	if (first_peak != values_end)
	{
		peaks.push_back(first_peak - values_begin);
		errors.push_back(get_half_value_width(
			first_peak - first_peak_index / 2,
			std::min(first_peak + first_peak_index / 2, values_end),
			first_peak
		));
	}

	auto prev_peak_it = first_peak;
	for (size_t i = 1; i < max_peaks; ++i)
	{
		auto const begin = std::min(prev_peak_it + (first_peak_index / 2), values_end);
		auto const end   = std::min(begin + first_peak_index, values_end);

		if (begin == values_end)
		{
			break;
		}

		auto peak_it = std::max_element(begin, end);
		peaks.push_back(peak_it - values_begin);
		prev_peak_it = peak_it;

		errors.push_back(get_half_value_width(begin, end, peak_it));
	}

	return result;
}

DLL_API void CDECL get_peaks(
	double in[], size_t in_size,
	double out_freq[], size_t out_freq_max_size,
	double out_hvw[], size_t out_hvw_max_size,
	double out_amp[], size_t out_amp_max_size,
	double df,
	double threshold,
	double min_search_freq, double max_search_freq
)
{
	auto const out_size = std::min(std::min(out_freq_max_size, out_hvw_max_size), out_amp_max_size);

	if (out_size == 0 || in_size == 0)
	{
		return;
	}

	std::vector<double> values;
	values.resize(in_size);

	std::transform(in, in + in_size, values.begin(), [](double val) {
		return 20 * std::log10(val);
	});

	auto const [peaks, errors] = find_peaks(values, out_size, df, threshold, min_search_freq, max_search_freq);

	for (size_t i = 0; i < out_size && i < peaks.size(); ++i)
	{
		out_freq[i] = peaks[i] * df;
		out_hvw[i]  = errors[i] * df;
		out_amp[i]  = in[peaks[i]];
	}
}

DLL_API void CDECL output_to_file(
	char path[],
	double freq[], size_t freq_size,
	double hvw[], size_t hvw_size,
	double amp[], size_t amp_size
)
{
	std::ofstream file(path, std::ios_base::app);
	if (!file.good())
	{
		return;
	}

	auto const size = std::min(std::min(freq_size, hvw_size), amp_size);

	file << "f	f_hvw	A\n";
	for (size_t i = 0; i < size; ++i)
	{
		file << freq[i] << '\t' << hvw[i] << '\t' << amp[i] << '\n';
	}
	file << '\n';
}

DLL_API uint32_t CDECL guess_letter(
	double in[], size_t in_size,
	double samples[], size_t samples_height, size_t samples_width
)
{
	std::vector<double> samples_norm(samples_height * samples_width, 0.0);

	for (size_t i = 0; i < samples_height; ++i)
	{
		auto const og_begin = samples + i * samples_width;
		auto const og_end   = og_begin + samples_width;

		auto const max_val = *std::max_element(og_begin, og_end);

		auto og_it = og_begin;
		auto it = samples_norm.begin() + i * samples_width;
		auto const end = it + samples_width;
		for (; it != end; ++it, ++og_it)
		{
			*it = *og_it / max_val;
		}
	}

	std::vector<double> in_norm(in_size, 0.0);

	{
		auto const in_begin = in;
		auto const in_end   = in + in_size;

		auto const max_val = *std::max_element(in_begin, in_end);

		auto in_it = in_begin;
		auto it = in_norm.begin();
		auto const end = in_norm.end();
		for (; it != end; ++it, ++in_it)
		{
			*it = *in_it / max_val;
		}
	}

	auto const get_sample = [&, samples_width, samples_height] (size_t i, size_t j)
	{
		if (i >= samples_height || j >= samples_width)
		{
			return 0.0;
		}
		return samples_norm[i * samples_width + j];
	};

	auto const calculate_error = [] (double x, double y)
	{
		auto const x_db = std::log10(x);
		auto const y_db = std::log10(y);
		return (x_db - y_db) * (x_db - y_db);
	};

	std::vector<double> errors(samples_height, 0.0);

	for (size_t i = 0; i < samples_height; ++i)
	{
		double error = 0.0;
		for (size_t j = 0; j < std::min(in_size, samples_width); ++j)
		{
			error += calculate_error(in_norm[j], get_sample(i, j));
		}
		errors[i] = error;
	}

	auto const min_it = std::min_element(errors.begin(), errors.end());
	return min_it - errors.begin();
}

#ifdef __GNUC__
#pragma GCC visibility pop
#endif // gcc
