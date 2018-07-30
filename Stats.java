import java.util.*;

/*
 * StatsJ Copyright 2018
 * Brian TerBush
 * brianterbush@gmail.com
 *
 * The StatsJ package is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (either version 2 of the License or, at your
 * option, any later version), provided that this notice and the name of the
 * author appear in all copies. Upon request to the author, some of the packages
 * in the Java Stats distribution can be licensed under the GNU Lesser General
 * Public License as published by the Free Software Foundation (either
 * version 2 of the License, or (at your option) any later version).
 * If you're using the software, please notify brianterbush@gmail.com so
 * that you can receive updates and patches. StatsJ is distributed
 * "as is", in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with the Java Stats distribution. If not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/**
 * The class <code>Stats</code> contains methods for performing basic
 * statistical operations that describe a sample or population, including
 * measures of central tendency or location, measures of variability or
 * dispersion, and measures of strength of relationships.
 * 
 * @author Brian TerBush
 */
public class Stats {

	protected Stats() {
	}

	/**
	 * Correlation is a measure of the degree to which the values contained in
	 * the arrays are related. If there is perfect relationship between the two
	 * arrays, we have a correlation of 1; if there is positive correlation,
	 * whenever one variable has a high (low) value, so does the other. If there
	 * is a perfect relationship with negative slope between the two variables,
	 * we have a correlation coefficient of -1; if there is negative
	 * correlation, whenever one variable has a high (low) value, the other has
	 * a low (high) value. A correlation coefficient of 0 means that there is no
	 * linear relationship between the variables.
	 * <p>
	 * The equation for correlation is:<br>
	 * <p>
	 * <img src="doc-files/cor.gif" alt="Correlation">
	 * <p>
	 * 
	 * @param xNumbers
	 *            the first array of numbers
	 * @param yNumbers
	 *            the second array of numbers
	 * @return measure between -1.0 and 1.0 of how strongly
	 *         <code>xNumbers</code> and <code>yNumbers</code> are related
	 * @exception ArithmeticException
	 *                thrown when <code>xNumbers</code> and
	 *                <code>yNumbers</code> are not the same length or the
	 *                standard deviation of <code>xNumbers</code> or
	 *                <code>yNumbers</code> is 0
	 * @exception NullPointerException
	 *                thrown if <code>xNumbers</code> or <code>yNumbers</code>
	 *                are null
	 * @see #cov(double[], double[])
	 */
	public static double corr(double[] xNumbers, double[] yNumbers)
			throws ArithmeticException {

		if (areEqual(xNumbers, yNumbers))
			return 1.0;

		double sx = popStdDev(xNumbers);

		if (sx == 0.0)
			throw new ArithmeticException(
					"The standard deviation of the values in xNumbers equals zero.");

		double sy = popStdDev(yNumbers);

		if (sy == 0.0)
			throw new ArithmeticException(
					"The standard deviation of the values in yNumbers equals zero.");

		return cov(xNumbers, yNumbers) / (sx * sy);
	}

	/**
	 * Covariance is the average of the products of deviations for each data
	 * point pair. Use covariance to determine the relationship between two data
	 * sets. A positive value means that the numbers are directly related. A
	 * negative value means that the variables are inversely related. A large
	 * magnitude indicates the likelihood of these associations are stronger.
	 * <p>
	 * The equation for covariance is:<br>
	 * <p>
	 * <img src="doc-files/cov.gif" alt="Covariance">
	 * <p>
	 * 
	 * @param xNumbers
	 *            the first array of numbers
	 * @param yNumbers
	 *            the second array of numbers
	 * @return the extent to which two arrays of numbers are related
	 * @exception ArithmeticException
	 *                thrown when <code>xNumbers</code> and
	 *                <code>yNumbers</code> are not the same length or either
	 *                is empty
	 * @exception NullPointerException
	 *                thrown when <code>xNumbers</code> or
	 *                <code>yNumbers</code> is null
	 */
	public static double cov(double[] xNumbers, double[] yNumbers)
			throws ArithmeticException {
		if (xNumbers.length == 0)
			throw new ArithmeticException(
					"ArrayNumbers must contain at least one value.");

		if (yNumbers.length == 0)
			throw new ArithmeticException(
					"Array yNumbers cannot have zero length.");

		if (xNumbers.length != yNumbers.length)
			throw new ArithmeticException(
					"Arrays xNumbers and y must have equal length.");

		double xMean = mean(xNumbers);
		double yMean = mean(yNumbers);
		double sum = 0.0;

		for (int i = 0; i < xNumbers.length; i++) {
			sum += (xNumbers[i] - xMean) * (yNumbers[i] - yMean);
		}

		return sum / (double) (xNumbers.length);
	}

	/**
	 * Frequency is a measure of the number of occurrences of
	 * <code>number</code> in <code>numbers</code>.
	 * <p>
	 * 
	 * @param numbers
	 *            numbers used to compute the frequency
	 * @param number
	 *            number to search for
	 * @return number of occurrences of a number in a data set
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static int freq(double[] numbers, double number) {
		int frequency = 0;

		for (int i = 0; i < numbers.length; i++) {
			if (numbers[i] == number)
				frequency++;
		}

		return frequency;
	}

	/**
	 * Returns the geometric mean of an array of positive data.
	 * <p>
	 * The equation for the geometric mean is:
	 * <p>
	 * <img src="doc-files/gmean.gif" alt="Geometric Mean">
	 * 
	 * 
	 * @param numbers
	 *            numbers used to compute the geometric mean
	 * @return sample geometric mean
	 * @exception ArithmeticException
	 *                thrown when any value in <code>y</code> is not positive
	 * @exception NullPointerException
	 *                thrown when <code>y</code> is null
	 */
	public static double geoMean(double[] numbers) throws ArithmeticException {
		if (numbers.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		double result = (numbers.length > 0) ? 1.0 : 0.0;

		for (int i = 0; i < numbers.length; i++) {
			if (numbers[i] <= 0.0)
				throw new ArithmeticException(
						"Array must contain only positive values.");
			result *= numbers[i];
		}

		return Math.pow(result, (1.0 / (double) numbers.length));
	}

	/**
	 * Frequency is a measure of the number of occurrences of
	 * <code>number</code> in <code>numbers</code>.
	 * <p>
	 * 
	 * @param numbers
	 *            numbers used to compute the frequency
	 * @param number
	 *            number to search for
	 * @return number of occurrences of a number in a data set
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	private static int freq(int[] numbers, int number) {
		int frequency = 0;

		for (int i = 0; i < numbers.length; i++) {
			if (numbers[i] == number)
				frequency++;
		}

		return frequency;
	}

	/**
	 * Rank represents the rank in descending order of a number in a list of
	 * numbers.
	 * <p>
	 * 
	 * @param numbers
	 *            numbers to analyze
	 * @param number
	 *            number to search for
	 * @return sorted rank of a number in a list of numbers
	 * @exception ArithmeticException
	 *                thrown if <code>numbers</code> does not contain
	 *                <code>number</code> or <code>numbers</code> has zero
	 *                length
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static int rank(double[] numbers, double number)
			throws ArithmeticException {
		return rank(numbers, number, true);
	}

	/**
	 * Rank represents the sorted rank of a number in a list of numbers.
	 * <p>
	 * 
	 * @param numbers
	 *            numbers to analyze
	 * @param number
	 *            number to search for
	 * @param useDescendingOrder
	 *            specifies whether to rank using descending or ascending order
	 * @return sorted rank of a number in a list of numbers
	 * @exception ArithmeticException
	 *                thrown if the number is not in <code>numbers</code> or
	 *                <code>numbers</code> has zero length
	 * @exception NullPointerException
	 *                thrown if <code>numbers</code> is null
	 */
	public static int rank(double[] numbers, double number,
			boolean useDescendingOrder) throws ArithmeticException {

		if (numbers.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		if (freq(numbers, number) == 0)
			throw new ArithmeticException("Number: " + number
					+ " not found in array.");

		double[] sorted = sort(numbers);
		int rank = 0;

		if (useDescendingOrder) {
			for (int i = sorted.length - 1; i >= 0; i--, rank++) {
				if (sorted[i] == number)
					break;
			}
		} else {
			for (int i = 0; i < sorted.length; i++, rank++) {
				if (sorted[i] == number)
					break;
			}
		}
		return rank;
	}

	/**
	 * Returns the skewness of a distribution. Skewness characterizes the degree
	 * of asymmetry of a distribution around its mean. Positive skewness
	 * indicates a distribution with an asymmetric tail extending toward more
	 * positive values. Negative skewness indicates a distribution with an
	 * asymmetric tail extending toward more negative values.
	 * <p>
	 * The equation for skewness is:
	 * <p>
	 * <img src="doc-files/skew.gif" alt="Skew">
	 * <p>
	 * 
	 * @param numbers
	 *            array of numbers used to compute the skewness
	 * @return the side with the tail (positive = right) and (negative = left)
	 * @exception ArithmeticException
	 *                thrown if <code>numbers</code> does not contain at least
	 *                three values or the sample standard deviation of
	 *                <code>numbers</code> is zero
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static double skew(double[] numbers) throws ArithmeticException {
		if (numbers.length < 3) {
			throw new ArithmeticException(
					"Array must contain at least three values.");
		}

		double s = stdDev(numbers);

		if (s == 0) {
			throw new ArithmeticException(
					"The sample standard deviation of the values in the array are zero.");
		}

		double u = mean(numbers);
		double sum = 0.0;

		for (int i = 0; i < numbers.length; i++) {
			sum += Math.pow((numbers[i] - u) / s, 3.0);
		}

		double n = (double) numbers.length;

		return n * sum / ((n - 1.0) * (n - 2.0));
	}

	/**
	 * Kurtosis is a measure of whether the data is peaked or flat relative to a
	 * normal distribution. That is, data sets with a high kurtosis tend to have
	 * a distinct peak near the mean, decline rather rapidly, and have heavy
	 * tails. Data sets with low kurtosis tend to have a flat top near the mean
	 * rather than a sharp peak.
	 * <p>
	 * The equation for the kurtosis is:
	 * <p>
	 * <img src="doc-files/kurt.gif" alt="Kurtosis">
	 * <p>
	 * 
	 * @param numbers
	 *            array of numbers used to compute the kurtosis
	 * @return measure of whether the values are peaked or flat relative to a
	 *         normal distribution
	 * @exception ArithmeticException
	 *                thrown when <code>numbers</code> has a length of zero
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static double kurt(double[] numbers) throws ArithmeticException {
		if (numbers.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		double u = mean(numbers);
		double s = stdDev(numbers);
		double sum = 0.0D;
		int n = numbers.length;

		if (numbers.length == 0 || s == 0.0D)
			return (0.0D / 0.0D);

		for (int i = 0; i < numbers.length; i++) {
			sum += Math.pow((numbers[i] - u) / s, 4.0);
		}

		return ((double) (n * (n + 1)) * sum)
				/ (double) ((n - 1) * (n - 2) * (n - 3))
				- (3.0 * Math.pow(n - 1, 2.0)) / (double) ((n - 2) * (n - 3));
	}

	/**
	 * Trimmed mean is a measure of variation excluding extreme values. A 5%
	 * trimmed mean ignores the largest and smallest 5% of values. This
	 * Eliminates extreme cases from the computation of the mean results in a
	 * better estimate of central tendency when the data is nonnormal.
	 * <p>
	 * 
	 * @param numbers
	 *            array of numbers used to compute the trimmed mean
	 * @return double returns a 5% trimmed mean
	 * @exception ArithmeticException
	 *                thrown if <code>numbers</code> has zero length
	 * @exception NullPointerException
	 *                thrown if <code>numbers</code> is null
	 */
	public static double trimmedMean(double[] numbers)
			throws ArithmeticException {
		if (numbers.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		double[] sorted = sort(numbers);
		double sum = 0.0;
		int start = (int) (sorted.length * 0.05);
		int stop = (int) (Math.ceil(sorted.length * 0.95));

		for (int i = start; i < stop; i++) {
			sum += sorted[i];
		}
		return sum / (stop - start);
	}

	/**
	 * Mean is a measure of location, commonly called the average. Its value
	 * depends equally on all of the data which may include outliers. It may not
	 * appear representative of the central region for skewed data sets. It is
	 * especially useful as being representative of the whole sample for use in
	 * subsequent calculations.
	 * <p>
	 * The equation for the mean is:
	 * <p>
	 * <img src="doc-files/mean.gif" alt="Mean">
	 * <p>
	 * 
	 * @param numbers
	 *            array of numbers used to compute the mean
	 * @return average value of the values contained in the array
	 * @exception ArithmeticException
	 *                thrown when <code>numbers</code> has a length of zero
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static double mean(double[] numbers) throws ArithmeticException {
		if (numbers.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		return sum(numbers) / (double) numbers.length;
	}

	/**
	 * Median is a measure of location that represents the value in the middle
	 * of the numbers after being sorted. It is generally a good descriptive
	 * measure of the location which works well for skewed data, or data with
	 * outliers.
	 * <p>
	 * 
	 * @param numbers
	 *            array of numbers used to compute the median
	 * @return value halfway through the array
	 * @exception ArithmeticException
	 *                thrown when <code>numbers</code> has a length of zero
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static double median(double[] numbers) throws ArithmeticException {
		if (numbers.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		double[] sorted = sort(numbers);

		int index = (int) (sorted.length / 2.0);

		// if odd
		if ((sorted.length % 2) == 1)
			return sorted[index];

		return (sorted[index - 1] + sorted[index]) / 2;
	}

	/**
	 * The mode of a set of observations is the specific value that occurs with
	 * the greatest frequency. There may be more than one mode in a set of
	 * observations, if there are several values that all occur with the
	 * greatest frequency. A mode may also not exist; this is true if all the
	 * observations occur with the same frequency.
	 * <p>
	 * 
	 * @param numbers
	 *            array of numbers used to compute the mode
	 * @return number the occurs most frequently
	 * @exception ArithmeticException
	 *                thrown if there is not a value that occurs most frequently
	 *                in the array or <code>numbers</code> has zero length
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static double[] mode(double[] numbers) throws ArithmeticException {
		if (numbers.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		double[] uniqueValues = uniqueValues(numbers);

		if (uniqueValues.length == 1)
			return new double[] { uniqueValues[0] };

		int[] counts = new int[uniqueValues.length];

		for (int i = 0; i < uniqueValues.length; i++) {
			counts[i] = freq(numbers, uniqueValues[i]);
		}

		int max = max(counts);
		int modeCount = freq(counts, max);

		if (modeCount == uniqueValues.length)
			return new double[] {};

		double[] modes = new double[modeCount];

		for (int i = 0, j = 0; i < counts.length; i++) {
			if (counts[i] == max) {
				modes[j++] = uniqueValues[i];
			}
		}

		return modes;
	}

	/**
	 * Range is a measure of difference between the largest and the smallest
	 * observed values in the array.
	 * <p>
	 * 
	 * @param numbers
	 *            array of numbers used to compute the range
	 * @return difference between the largest and smallest values in the array
	 * @exception ArithmeticException
	 *                thrown when <code>numbers</code> has a length less than
	 *                zero.
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static double range(double[] numbers) throws ArithmeticException {
		if (numbers.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		if (numbers.length == 1)
			return 0.0;

		return max(numbers) - min(numbers);
	}

	/**
	 * Estimates standard deviation based on a sample. The standard deviation is
	 * a measure of how widely values are dispersed from the mean.
	 * <p>
	 * The equation for variance is:
	 * <p>
	 * <img src="doc-files/sdev.gif" alt="Standard Deviation">
	 * <p>
	 * 
	 * @param numbers
	 *            array of numbers used to compute the standard deviation
	 * @return dispersion of a set of data
	 * @exception ArithmeticException
	 *                thrown when <code>numbers</code> has a length of zero
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static double stdDev(double[] numbers) throws ArithmeticException {
		return Math.sqrt(var(numbers));
	}

	/**
	 * Calculates standard deviation based on the entire population given as
	 * arguments. The standard deviation is a measure of how widely values are
	 * dispersed from the mean.
	 * <p>
	 * The equation for variance is:
	 * <p>
	 * <img src="sdevp.gif" alt="Population Standard Deviation">
	 * <p>
	 * 
	 * @param numbers
	 *            array of numbers used to compute the standard deviation
	 * @return measure of how widely values are dispersed from the mean
	 * @exception ArithmeticException
	 *                thrown when <code>numbers</code> has a length of zero
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static double popStdDev(double[] numbers) throws ArithmeticException {
		return Math.sqrt(popVar(numbers));
	}

	/**
	 * Computes the sum of all the values in the array.
	 * 
	 * @param numbers to sum
	 * @return sum of the values
	 */
	public static double sum(double[] numbers) {
		double sum = 0;

		for (int i = 0; i < numbers.length; i++) {
			sum += numbers[i];
		}

		return sum;
	}

	/**
	 * Computes the sum of the values in the array.
	 * 
	 * @return sum of squares
	 */
	private static double sumOfSquares(double[] numbers) {
		double sumsqr = 0;
		for (int i = 0; i < numbers.length; i++) {
			sumsqr += numbers[i] * numbers[i];
		}
		return sumsqr;
	}

	/**
	 * Sample variance is a measure of the spread of or dispersion of the values
	 * in a sample population.
	 * <p>
	 * The equation for variance is:
	 * <p>
	 * <img src="doc-files/var.gif" alt="Variance">
	 * <p>
	 * 
	 * @param numbers
	 *            samples of the population
	 * @return estimate of variance based on a sample population
	 * @exception ArithmeticException
	 *                thrown when <code>numbers</code> has a length of zero
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static double var(double[] numbers) throws ArithmeticException {
		if (numbers.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		double sum = sum(numbers);

		return (sumOfSquares(numbers) - (sum * sum / numbers.length))
				/ (numbers.length - 1);
	}

	/**
	 * Population variance is a measure of the spread of or dispersion of the
	 * values in entire population.
	 * <p>
	 * The equation for population variance is:
	 * <p>
	 * <img src="doc-files/varp.gif" alt="Population Variance">
	 * <p>
	 * 
	 * @param numbers
	 *            samples of the population
	 * @return estimate of variance based on the entire population
	 */
	public static double popVar(double[] numbers) throws ArithmeticException {
		if (numbers.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		double sum = sum(numbers);

		return (sumOfSquares(numbers) - (sum * sum / numbers.length))
				/ (numbers.length);
	}

	/**
	 * Returns the smallest number in the array.
	 * <p>
	 * 
	 * @param numbers
	 *            array of numbers used to compute the min
	 * @return largest number found
	 * @exception ArithmeticException
	 *                thrown when <code>numbers</code> has a length of zero
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static double min(double[] numbers) throws ArithmeticException {
		if (numbers.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		double min = Double.MAX_VALUE;

		for (int i = 0; i < numbers.length; i++) {
			if (numbers[i] == Double.MIN_VALUE)
				return Double.MIN_VALUE;

			if (numbers[i] < min)
				min = numbers[i];
		}

		return min;
	}

	/**
	 * Returns the largest number in the array.
	 * <p>
	 * 
	 * @param numbers
	 *            array of numbers used to compute the max
	 * @return largest number found
	 * @exception ArithmeticException
	 *                thrown when <code>numbers</code> has a length of zero
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	public static double max(double[] numbers) throws ArithmeticException {
		if (numbers.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		double max = Double.MIN_VALUE;

		for (int i = 0; i < numbers.length; i++) {
			if (numbers[i] == Double.MAX_VALUE)
				return Double.MAX_VALUE;

			if (numbers[i] > max)
				max = numbers[i];
		}

		return max;
	}

	/**
	 * Returns the largest number in the array.
	 * <p>
	 * 
	 * @param numbers
	 *            array of numbers used to compute the max
	 * @return largest number found
	 * @exception ArithmeticException
	 *                thrown when <code>numbers</code> has a length of zero
	 * @exception NullPointerException
	 *                thrown when <code>numbers</code> is null
	 */
	private static int max(int[] x) throws ArithmeticException {
		if (x.length == 0)
			throw new ArithmeticException(
					"Array must contain at least one value.");

		int max = Integer.MIN_VALUE;

		for (int i = 0; i < x.length; i++) {
			if (x[i] == Integer.MAX_VALUE)
				return Integer.MAX_VALUE;

			if (x[i] > max)
				max = x[i];
		}

		return max;
	}

	private static double[] uniqueValues(double[] arr) {
		Set set = new HashSet();

		for (int i = 0; i < arr.length; i++) {
			set.add(new Double(arr[i]));
		}

		double[] values = new double[set.size()];
		Iterator itr = set.iterator();
		int index = 0;

		while (itr.hasNext()) {
			values[index++] = ((Double) itr.next()).doubleValue();
		}

		return values;
	}

	private static boolean areEqual(double[] x, double[] y)
			throws ArithmeticException {
		if (x.length != y.length)
			throw new ArithmeticException("arrays must be same length.");

		for (int i = 0; i < x.length; i++) {
			if (x[i] != y[i])
				return false;
		}

		return true;
	}

	protected static double[] sort(double[] arr) {
		double[] sorted = new double[arr.length];
		System.arraycopy(arr, 0, sorted, 0, arr.length);

		if (sorted.length > 1)
			quickSort(sorted, 0, sorted.length - 1);

		return sorted;
	}

	protected static void quickSort(double arr[], int p, int r) {
		int k = p;
		int l = r;
		if (k < l) {
			double val = arr[(k + l) / 2];
			while (k <= l) {
				while (l > p && arr[l] > val)
					l--;
				while (k < r && arr[k] < val)
					k++;
				if (k <= l) {
					double d = arr[k];
					arr[k++] = arr[l];
					arr[l--] = d;
				}
			}
			if (p < l)
				quickSort(arr, p, l);
			if (k < r)
				quickSort(arr, k, r);
		}
	}
}