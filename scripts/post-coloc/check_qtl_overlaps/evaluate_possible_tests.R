require(metap)
require(sensitivitymv)
require(qvalue)
require(plotly)
set.seed(0)

# Assume a million GWAS hits,
# each tested for eQTLs in num_tissues tissues.

# Training set: will be used for choosing parameters on the more difficult methods
# Test set: will be used for evaluating performance

#num_trials = 1000000
num_trials = 10000
num_tissues = 40

power = array(0, dim=c(10,10,8))
significance = array(0, dim=c(10,10,8))
for (num_hits in 1:10)
{
	for (ss in 1:10)
	{
		signal_strength = (-10:-1)[ss]
		print("Starting next set:")
		print(paste(num_hits, signal_strength))

		# Negative set: a bunch of p-values distributed uniformly at random
		test_null = array(0, dim=c(num_trials, num_tissues))
		for (i in 1:num_tissues)
		{
			test_null[1:num_trials, i] = runif(num_trials) 
		}

		# Positive set: throw in a few actual signals
		test_positive = array(0, dim=c(num_trials, num_tissues))
		for (i in 1:num_tissues)
		{
			test_positive[1:num_trials, i] = runif(num_trials) 
			test_positive[,1:num_hits] = 10^(runif(num_trials*num_hits, signal_strength, -1))
		}

		# Now we'll test every method for power and for false positive rate trying to achieve p < 0.05

		# Method 1 (baseline): Select p-value uniformly at random.
		# This is one of the worst methods we could imagine, so 
		# any good method will perform better.

		test_null_results = runif(num_trials)
		test_positive_results = runif(num_trials)
		print("Baseline")
		print(sum(test_null_results < 0.05) / num_trials)		# False positive rate (should be around 0.05 if we trust the test)
		print(sum(test_positive_results < 0.05) / num_trials)		# Power to reject the null
		significance[num_hits, ss, 1] = sum(test_null_results < 0.05) / num_trials
		power[num_hits, ss, 1] = sum(test_positive_results < 0.05) / num_trials

		# Method 2 (min p-val): Choose the minimum p-value.
		# This invalid method will have great power, if we
		# don't mind having an enormous false positive rate.
		test_null_results = apply(test_null, 1, min)
		test_positive_results = apply(test_positive, 1, min)
		print("Min p-val")
		print(sum(test_null_results < 0.05) / num_trials)		# False positive rate (should be around 0.05 if we trust the test)
		print(sum(test_positive_results < 0.05) / num_trials)		# Power to reject the null
		significance[num_hits, ss, 2] = sum(test_null_results < 0.05) / num_trials
		power[num_hits, ss, 2] = sum(test_positive_results < 0.05) / num_trials

		# Method 3 (Bonferroni): Choose the minimum p-value.
		test_null_results = apply(test_null, 1, min) * num_tissues
		test_positive_results = apply(test_positive, 1, min) * num_tissues
		print("Bonferroni")
		print(sum(test_null_results < 0.05) / num_trials)		# False positive rate (should be around 0.05 if we trust the test)
		print(sum(test_positive_results < 0.05) / num_trials)		# Power to reject the null
		significance[num_hits, ss, 3] = sum(test_null_results < 0.05) / num_trials
		power[num_hits, ss, 3] = sum(test_positive_results < 0.05) / num_trials

		# Method 4 (BH): Correct with BH, then choose minimum p-value
		test_null_results = apply(test_null, 1, function(x)
			{
			       min(p.adjust(x, method="BH"))
			})
		test_positive_results = apply(test_positive, 1, function(x)
			{
			       min(p.adjust(x, method="BH"))
			})
		print("Benjamini-Hochberg")
		print(sum(test_null_results < 0.05) / num_trials)		# False positive rate (should be around 0.05 if we trust the test)
		print(sum(test_positive_results < 0.05) / num_trials)		# Power to reject the null
		significance[num_hits, ss, 4] = sum(test_null_results < 0.05) / num_trials
		power[num_hits, ss, 4] = sum(test_positive_results < 0.05) / num_trials

		# Method 5 (Binomial test)
		# Will require parameter tuning: creep down until we get a good value here
		cutpoint = 0.1
		test_null_results = apply(test_null, 1, function(x)
			{
				n = sum(x < cutpoint)
				1-pbinom(n, num_tissues, cutpoint)
			})
		test_positive_results = apply(test_positive, 1, function(x)
			{
				n = sum(x < cutpoint)
				1-pbinom(n, num_tissues, cutpoint)
			})
		print("Binomial")
		print(sum(test_null_results < 0.05) / num_trials)		# False positive rate (should be around 0.05 if we trust the test)
		print(sum(test_positive_results < 0.05) / num_trials)		# Power to reject the null
		significance[num_hits, ss, 5] = sum(test_null_results < 0.05) / num_trials
		power[num_hits, ss, 5] = sum(test_positive_results < 0.05) / num_trials

		# Method 6 (Fisher's combined)
		test_null_results = apply(test_null, 1, function(x)
			{
				return(sumlog(x)$p)
			})
		test_positive_results = apply(test_positive, 1, function(x)
			{
				return(sumlog(x)$p)
			})
		print("Fisher's combined")
		print(sum(test_null_results < 0.05) / num_trials)		# False positive rate (should be around 0.05 if we trust the test)
		print(sum(test_positive_results < 0.05) / num_trials)		# Power to reject the null
		significance[num_hits, ss, 6] = sum(test_null_results < 0.05) / num_trials
		power[num_hits, ss, 6] = sum(test_positive_results < 0.05) / num_trials

		# Method 7 (Zaykin's truncated product)
		# Will require parameter tuning
		test_null_results = apply(test_null, 1, function(x)
			{
				return(truncatedP(x))
			})
		test_positive_results = apply(test_positive, 1, function(x)
			{
				return(truncatedP(x))
			})
		print("Zaykin's truncated")
		print(sum(test_null_results < 0.05) / num_trials)		# False positive rate (should be around 0.05 if we trust the test)
		print(sum(test_positive_results < 0.05) / num_trials)		# Power to reject the null
		significance[num_hits, ss, 7] = sum(test_null_results < 0.05) / num_trials
		power[num_hits, ss, 7] = sum(test_positive_results < 0.05) / num_trials

		# Method 8 (KS test)
		test_null_results = apply(test_null, 1, function(x)
			{
				return(ks.test(x,punif)$p.value)
			})
		test_positive_results = apply(test_positive, 1, function(x)
			{
				return(ks.test(x,punif)$p.value)
			})
		print("Kolmogorov-Smirnov against uniform")
		print(sum(test_null_results < 0.05) / num_trials)		# False positive rate (should be around 0.05 if we trust the test)
		print(sum(test_positive_results < 0.05) / num_trials)		# Power to reject the null
		significance[num_hits, ss, 8] = sum(test_null_results < 0.05) / num_trials
		power[num_hits, ss, 8] = sum(test_positive_results < 0.05) / num_trials
	}
}

methods = c("Random", "Minimum p-value", "Bonferroni", "Benjamini-Hochberg", "Binomial (#values < 0.1)", "Fisher's combined", "Zaykin's truncated (trunc = 0.2)", "Kolmogorov-Smirnov against uniform")
for (method in 1:8)
{
	p = plot_ly(x=-10:-1, y=1:5, z=power[,,method], type="contour",
		      contours = list(
				          start = 0,
					      end = 1,
					      size = 0.05
					        )) %>% 
			layout(
		           	title = methods[method],
			       	xaxis = list(title = "log10(min P-value)"), yaxis = list(title = "number of true effects (out of 40)"))
	orca(p, paste0('images/method', method, 'contour.png'))

	p = plot_ly(x=-10:-1, y=1:10, z=significance[,,method], type="contour",
		      contours = list(
				          start = 0,
					      end = 1,
					      size = 0.05
					        )) %>% 
			layout(
		           	title = methods[method],
			       	xaxis = list(title = "log10(min P-value)"), yaxis = list(title = "number of true effects (out of 40)"))
	orca(p, paste0('images/method', method, 'sig_contour.png'))

}

