CLEANFILES = \
	*~ \
	markov-chain-sampler_TEST.hdf5 \
	markov-chain-sampler_TEST_density.hdf5 \
	pmc_sampler_TEST-mcmc-prerun.hdf5 \
	pmc_sampler_TEST-density.hdf5 \
	pmc_sampler_TEST-density-prerun.hdf5 \
	pmc_sampler_TEST-output.hdf5 \
	pmc_sampler_TEST-output-components.hdf5 \
	pmc_sampler_TEST-output-hc.hdf5 \
	pmc_sampler_TEST-output-resume.hdf5 \
	pmc_sampler_TEST-output-split.hdf5 \
	prior-sampler_TEST.hdf5 \
	proposal-functions_TEST-rdwr.hdf5 \
	proposal-functions_TEST-block-decomposition.hdf5
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

lib_LTLIBRARIES = libeosstatistics.la
libeosstatistics_la_SOURCES = \
	analysis.cc analysis.hh analysis-fwd.hh \
	chain-group.cc chain-group.hh \
	chi-squared.hh chi-squared.cc \
	density-wrapper.cc density-wrapper.hh \
	hierarchical-clustering.cc hierarchical-clustering.hh \
	histogram.cc histogram.hh \
	log-likelihood.cc log-likelihood.hh log-likelihood-fwd.hh \
	log-prior.cc log-prior.hh log-prior-fwd.hh \
	markov-chain.cc markov-chain.hh \
	markov-chain-sampler.cc markov-chain-sampler.hh \
	prior-sampler.cc prior-sampler.hh \
	proposal-functions.cc proposal-functions.hh \
	rvalue.cc rvalue.hh \
	simple-parameters.cc simple-parameters.hh \
	test-statistic.cc test-statistic.hh test-statistic-impl.hh \
	welford.cc welford.hh
libeosstatistics_la_LIBADD = -lpthread -lgsl -lgslcblas -lm -lMinuit2

if EOS_ENABLE_PMC
libeosstatistics_la_SOURCES += \
	population-monte-carlo-sampler.cc population-monte-carlo-sampler.hh
libeosstatistics_la_LIBADD += -lpmc -ldl
endif

include_eos_statisticsdir = $(includedir)/eos/statistics
include_eos_statistics_HEADERS = \
	analysis.hh analysis-fwd.hh \
	chain-group.hh \
	chi-squared.hh \
	density-wrapper.hh \
	hierarchical-clustering.hh \
	histogram.hh \
	log-likelihood.hh log-likelihood-fwd.hh \
	log-prior.hh log-prior-fwd.hh \
	markov-chain.hh \
	markov-chain-sampler.hh \
	prior-sampler.hh \
	proposal-functions.hh \
	rvalue.hh \
	simple-parameters.hh \
	welford.hh

if EOS_ENABLE_PMC
include_eos_statistics_HEADERS += population-monte-carlo-sampler.hh
endif

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	analysis_TEST \
	chi-squared_TEST \
	density-wrapper_TEST \
	hierarchical-clustering_TEST \
	histogram_TEST \
	log-likelihood_TEST \
	log-prior_TEST \
	markov-chain_TEST \
	markov-chain-sampler_TEST \
	prior-sampler_TEST \
	proposal-functions_TEST \
	rvalue_TEST \
	simple-parameters_TEST \
	welford_TEST
LDADD = \
	-lhdf5 \
	-lMinuit2 \
	$(top_builddir)/test/libeostest.a \
	libeosstatistics.la \
	$(top_builddir)/eos/libeos.la \
	$(top_builddir)/eos/utils/libeosutils.la

if EOS_ENABLE_PMC
TESTS += population-monte-carlo-sampler_TEST
LDADD += -lpmc -ldl
endif

check_PROGRAMS = $(TESTS)

analysis_TEST_SOURCES = analysis_TEST.cc analysis_TEST.hh

chi_squared_TEST_SOURCES = chi-squared_TEST.cc

density_wrapper_TEST_SOURCES = density-wrapper_TEST.cc density-wrapper_TEST.hh

hierarchical_clustering_TEST_SOURCES = hierarchical-clustering_TEST.cc

histogram_TEST_SOURCES = histogram_TEST.cc

log_likelihood_TEST_SOURCES = log-likelihood_TEST.cc

log_prior_TEST_SOURCES = log-prior_TEST.cc

markov_chain_TEST_SOURCES = markov-chain_TEST.cc

markov_chain_sampler_TEST_SOURCES = markov-chain-sampler_TEST.cc density-wrapper_TEST.cc

if EOS_ENABLE_PMC
population_monte_carlo_sampler_TEST_SOURCES = population-monte-carlo-sampler_TEST.cc density-wrapper_TEST.cc
endif

prior_sampler_TEST_SOURCES = prior-sampler_TEST.cc

proposal_functions_TEST_SOURCES = proposal-functions_TEST.cc

rvalue_TEST_SOURCES = rvalue_TEST.cc

simple_parameters_TEST_SOURCES = simple-parameters_TEST.cc

welford_TEST_SOURCES = welford_TEST.cc
