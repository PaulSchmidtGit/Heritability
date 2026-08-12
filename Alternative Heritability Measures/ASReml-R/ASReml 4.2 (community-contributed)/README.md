## ASReml-R 4.2.x versions (community-contributed)

The scripts in the **parent folder** were written for **ASReml-R 3.0**. The
scripts in *this* folder are the same examples adapted to the newer
**ASReml-R 4.2.x** syntax.

> ⚠️ **These updated scripts were not authored or tested by us.**
> They were kindly contributed by GitHub user
> [@mpeddu](https://github.com/mpeddu) in
> [pull request #5](https://github.com/PaulSchmidtGit/Heritability/pull/5).
> The repository maintainer does not have an ASReml licence and therefore
> **could not run or verify this code**. Please review and test it yourself
> before relying on the results.

The main changes from ASReml-R 3.0 to 4.2.x reflected here are:

* `predict()` no longer nests its output under `$pred` (e.g. `predict(...)$pred$sed` → `predict(...)$sed`),
* variance-component labels changed (e.g. `varcomp['gen!gen.var', ]` → `varcomp['gen', ]`),
* `predict()` returns S4 matrices, so `reshape::melt()` needs an explicit `as.matrix()`.

---

📦 **Looking for a ready-made R package?** Most of these heritability measures
are now available in the **[`heritable`](https://github.com/anu-aagi/heritable)**
R package (docs: <https://anu-aagi.github.io/heritable/>).
