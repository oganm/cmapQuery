
Readme
------

This is currently more of a "notes to self" rather than a readme.

This repository aims to replicate the analysis performed at Lamb et al 2006. I but allow using the L1000 data instead of the original Cmap data. With L1000, the data is much larger and there are a few different version of the stage 5 data around created by different approaches

-   The real stage 5 data by LINCS people
-   Our stage 5 data calculated by comparing to controls
-   Avi Mayan's stage 5 data, calculated by comparing to controls from the same lane. it's a little confusing. comunications pending

Since I want to experiment with all of those, this needs to be dataset independent. However my original methodology relied on memoisation for speed which is problematic if the background data can change. I need to decide when to do this memoisation and how much should I expose in on user side

Current process is to run `preCalcRandomKs` before running the experiment manually to get a memoised function and use it as an input to `connectivityMapEnrichment`. This will take a long time for the new L1000 data. I can start to see why Mayan's group abondoned this method.

To be able to quickly compare the studies it might make sense to use Mayan's web apps. I am implementing wrappers for the API calls of the apps. Currently the FWD site is not functional from their side
