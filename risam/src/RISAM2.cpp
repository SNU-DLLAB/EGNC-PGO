/* ----------------------------------------------------------------------------
 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

#include "risam/RISAM2.h"

#include <gtsam/nonlinear/ISAM2-impl.h>

#include <boost/math/distributions/chi_squared.hpp>

#include "risam/GraduatedFactor.h"

namespace risam {

/**
 * ####  ######     ###    ##     ##  #######
 *  ##  ##    ##   ## ##   ###   ### ##     ##
 *  ##  ##        ##   ##  #### ####        ##
 *  ##   ######  ##     ## ## ### ##  #######
 *  ##        ## ######### ##     ## ##
 *  ##  ##    ## ##     ## ##     ## ##
 * ####  ######  ##     ## ##     ## #########
 */

/* ************************************************************************* */
void RISAM2::set_mode(int mode) {
  switch (mode) {
    case 0:
      updateMuInitOnConvergence = std::bind(&RISAM2::updateMuInitOnConvergence_m0, this);
      relinearizeFactor = std::bind(&RISAM2::relinearizeFactor_m0, this, std::placeholders::_1, std::placeholders::_2,
                                    std::placeholders::_3);
      pushBackFactors = std::bind(&RISAM2::pushBackFactors_m0, this, std::placeholders::_1, std::placeholders::_2,
                                  std::placeholders::_3, std::placeholders::_4, std::placeholders::_5,
                                  std::placeholders::_6, std::placeholders::_7, std::placeholders::_8);
      break;
    case 1:
      updateMuInitOnConvergence = std::bind(&RISAM2::updateMuInitOnConvergence_m1, this);
      relinearizeFactor = std::bind(&RISAM2::relinearizeFactor_m1, this, std::placeholders::_1, std::placeholders::_2,
                                    std::placeholders::_3);
      pushBackFactors = std::bind(&RISAM2::pushBackFactors_m1, this, std::placeholders::_1, std::placeholders::_2,
                                  std::placeholders::_3, std::placeholders::_4, std::placeholders::_5,
                                  std::placeholders::_6, std::placeholders::_7, std::placeholders::_8);
      break;
    case 2:
      updateMuInitOnConvergence = std::bind(&RISAM2::updateMuInitOnConvergence_m2, this);
      relinearizeFactor = std::bind(&RISAM2::relinearizeFactor_m2, this, std::placeholders::_1, std::placeholders::_2,
                                    std::placeholders::_3);
      pushBackFactors = std::bind(&RISAM2::pushBackFactors_m2, this, std::placeholders::_1, std::placeholders::_2,
                                  std::placeholders::_3, std::placeholders::_4, std::placeholders::_5,
                                  std::placeholders::_6, std::placeholders::_7, std::placeholders::_8);
      break;
    case 3:
      updateMuInitOnConvergence = std::bind(&RISAM2::updateMuInitOnConvergence_m3, this);
      relinearizeFactor = std::bind(&RISAM2::relinearizeFactor_m3, this, std::placeholders::_1, std::placeholders::_2,
                                    std::placeholders::_3);
      pushBackFactors = std::bind(&RISAM2::pushBackFactors_m3, this, std::placeholders::_1, std::placeholders::_2,
                                  std::placeholders::_3, std::placeholders::_4, std::placeholders::_5,
                                  std::placeholders::_6, std::placeholders::_7, std::placeholders::_8);
      break;
    default:
      throw std::invalid_argument("Invalid mode");
  }
}

/* ************************************************************************* */
ISAM2Result RISAM2::update(const NonlinearFactorGraph& newFactors, const Values& newTheta, bool known_inliers,
                           const FactorIndices& removeFactorIndices,
                           const boost::optional<FastMap<Key, int>>& constrainedKeys,
                           const boost::optional<FastList<Key>>& noRelinKeys,
                           const boost::optional<FastList<Key>>& extraReelimKeys, bool force_relinearize) {
  ISAM2UpdateParams params;
  params.constrainedKeys = constrainedKeys;
  params.extraReelimKeys = extraReelimKeys;
  params.force_relinearize = force_relinearize;
  params.noRelinKeys = noRelinKeys;
  params.removeFactorIndices = removeFactorIndices;

  last_update_info_ = RiSAMUpdateInfo();

  new_factors_include_gnc_ = false;
  for (auto& factor : newFactors) {
    if (boost::dynamic_pointer_cast<GraduatedFactor>(factor)) {
      new_factors_include_gnc_ = true;
    }
  }

  auto result = update(newFactors, newTheta, known_inliers, params);

  if (risam_params_.converge_after_new_gnc && new_factors_include_gnc_) {
    iterateToConvergence();
  }

  // If we have converged update mu_inits_based on status
  if (risam_params_.increment_outlier_mu && convex_factors_.size() == 0 &&
      delta_.norm() / delta_.size() < risam_params_.outlier_mu_avg_var_convergence_thresh) {
    updateMuInitOnConvergence();
  }

  return result;
}

/* ************************************************************************* */
ISAM2Result RISAM2::update(const NonlinearFactorGraph& newFactors, const Values& newTheta, bool known_inliers,
                           const ISAM2UpdateParams& updateParams) {
  this->update_count_ += 1;
  UpdateImpl::LogStartingUpdate(newFactors, *this);
  ISAM2Result result(params_.enableDetailedResults);
  UpdateImpl update(params_, updateParams);

  // 1. Add any new factors \Factors:=\Factors\cup\Factors'.
  pushBackFactors(newFactors, known_inliers, &nonlinearFactors_, &linearFactors_, &variableIndex_,
                  &result.newFactorsIndices, &result.keysWithRemovedFactors, updateParams);
  update.computeUnusedKeys(newFactors, variableIndex_, result.keysWithRemovedFactors, &result.unusedKeys);

  // 2. Initialize any new variables \Theta_{new} and add
  // \Theta:=\Theta\cup\Theta_{new}.
  addVariables(newTheta, result.details());
  for (auto& key : newTheta.keys()) variable_ordering_.push_back(key);
  if (params_.evaluateNonlinearError) update.error(nonlinearFactors_, calculateEstimate(), &result.errorBefore);

  // 3. Mark linear update
  update.gatherInvolvedKeys(newFactors, nonlinearFactors_, result.keysWithRemovedFactors, &result.markedKeys);
  KeySet updateDirectlyInvolvedKeys = result.markedKeys;
  update.updateKeys(result.markedKeys, &result);

  KeySet relinKeys;
  result.variablesRelinearized = 0;
  if (update.relinarizationNeeded(update_count_)) {
    // 4. Mark keys in \Delta above threshold \beta:
    relinKeys = update.gatherRelinearizeKeys(roots_, delta_, fixedVariables_, &result.markedKeys);
    update.recordRelinearizeDetail(relinKeys, result.details());
    if (!relinKeys.empty()) {
      // 5. Mark cliques that involve marked variables \Theta_{J} and ancestors.
      update.findFluid(roots_, relinKeys, &result.markedKeys, result.details());
      // 6. Update linearization point for marked variables:
      // \Theta_{J}:=\Theta_{J}+\Delta_{J}.
      UpdateImpl::ExpmapMasked(delta_, relinKeys, &theta_);
    }
    result.variablesRelinearized = result.markedKeys.size();
  }

  // 7. Linearize new factors
  update.linearizeNewFactors(newFactors, theta_, nonlinearFactors_.size(), result.newFactorsIndices, &linearFactors_);
  update.augmentVariableIndex(newFactors, result.newFactorsIndices, &variableIndex_);

  // 8. Redo top of Bayes tree and update data structures
  recalculate(updateParams, relinKeys, &result, updateDirectlyInvolvedKeys);
  if (!result.unusedKeys.empty()) removeVariables(result.unusedKeys);
  result.cliques = this->nodes().size();

  if (params_.evaluateNonlinearError) update.error(nonlinearFactors_, calculateEstimate(), &result.errorAfter);

  // 9. Update Delta, this will force us to do a back sub.
  // Technically we could skip this, but usually people call getEstimate after each update anyways
  // If you do that it actually is slower to avoid calling delta because it gets called > #updates
  updateDelta();
  return result;
}

/* ************************************************************************* */
void RISAM2::recalculate(const ISAM2UpdateParams& updateParams, const KeySet& relinKeys, ISAM2Result* result,
                         const KeySet& updateDirectlyInvolvedKeys) {
  UpdateImpl::LogRecalculateKeys(*result);

  /** (riSAM)
   * Now we need to actually preform the bayes tree modification including keys marked as convex.
   * The set of all convex factors is included in class field convex_factors_.
   */
  auto convexKeys = KeyVector();
  for (auto factor_idx : convex_factors_) {
    convexKeys.insert(convexKeys.end(), nonlinearFactors_[factor_idx]->begin(), nonlinearFactors_[factor_idx]->end());
  }

  if (!result->markedKeys.empty() || !result->observedKeys.empty() || !convexKeys.empty()) {
    /** (riSAM)
     * We need to track specially all the keys that get affected by the new factors so that we can convexify these
     * factors. To do so we need to traverse the tree so we get all factors along path to root.
     */
    auto updateInvolvedKeys =
        traverseTop(KeyVector(updateDirectlyInvolvedKeys.begin(), updateDirectlyInvolvedKeys.end()));

    KeyVector involvedKeys;
    involvedKeys.insert(involvedKeys.end(), result->markedKeys.begin(), result->markedKeys.end());
    involvedKeys.insert(involvedKeys.end(), convexKeys.begin(), convexKeys.end());

    // Remove top of Bayes tree and convert to a factor graph:
    // (a) For each affected variable, remove the corresponding clique and all
    // parents up to the root. (b) Store orphaned sub-trees \BayesTree_{O} of
    // removed cliques.
    GaussianBayesNet affectedBayesNet;
    Cliques orphans;
    this->removeTop(involvedKeys, &affectedBayesNet, &orphans);

    // FactorGraph<GaussianFactor> factors(affectedBayesNet);
    // bug was here: we cannot reuse the original factors, because then the
    // cached factors get messed up [all the necessary data is actually
    // contained in the affectedBayesNet, including what was passed in from the
    // boundaries, so this would be correct; however, in the process we also
    // generate new cached_ entries that will be wrong (ie. they don't contain
    // what would be passed up at a certain point if batch elimination was done,
    // but that's what we need); we could choose not to update cached_ from
    // here, but then the new information (and potentially different variable
    // ordering) is not reflected in the cached_ values which again will be
    // wrong] so instead we have to retrieve the original linearized factors AND
    // add the cached factors from the boundary

    // ordering provides all keys in conditionals, there cannot be others
    // because path to root included
    FastList<Key> affectedKeys;
    for (const auto& conditional : affectedBayesNet)
      affectedKeys.insert(affectedKeys.end(), conditional->beginFrontals(), conditional->endFrontals());

    // (riSAM Debug)
    if (RISAM_DEBUG) {
      std::cout << "Update: [";
      for (auto k : updateInvolvedKeys) std::cout << DefaultKeyFormatter(k) << ", ";
      std::cout << "]" << std::endl;

      std::cout << "Convex: [";
      for (auto k : affectedKeys) std::cout << DefaultKeyFormatter(k) << ", ";
      std::cout << "]" << std::endl;
    }
    last_update_info_.updateInvolvedKeys = updateInvolvedKeys;
    last_update_info_.affectedKeysConvex.insert(affectedKeys.begin(), affectedKeys.end());

    KeySet affectedKeysSet;
    static const double kBatchThreshold = 0.65;
    if (affectedKeys.size() >= theta_.size() * kBatchThreshold) {
      // Do a batch step - reorder and relinearize all variables
      recalculateBatch(updateParams, updateInvolvedKeys, &affectedKeysSet, result);
    } else {
      recalculateIncremental(updateParams, relinKeys, affectedKeys, updateInvolvedKeys, &affectedKeysSet, &orphans,
                             result);
    }

    // Root clique variables for detailed results
    if (result->detail && params_.enableDetailedResults) {
      for (const auto& root : roots_)
        for (Key var : *root->conditional()) result->detail->variableStatus[var].inRootClique = true;
    }

    // Update replaced keys mask (accumulates until back-substitution happens)
    deltaReplacedMask_.insert(affectedKeysSet.begin(), affectedKeysSet.end());
  }
}

/* ************************************************************************* */
void RISAM2::recalculateBatch(const ISAM2UpdateParams& updateParams, const KeySet& updateInvolvedKeys,
                              KeySet* affectedKeysSet, ISAM2Result* result) {
  br::copy(variableIndex_ | br::map_keys, std::inserter(*affectedKeysSet, affectedKeysSet->end()));

  // Removed unused keys:
  VariableIndex affectedFactorsVarIndex = variableIndex_;

  affectedFactorsVarIndex.removeUnusedVariables(result->unusedKeys.begin(), result->unusedKeys.end());

  for (const Key key : result->unusedKeys) {
    affectedKeysSet->erase(key);
  }

  Ordering order;
  if (risam_params_.use_custom_ordering_constraint) {
    FastMap<Key, int> constraintGroups;
    if (updateInvolvedKeys.size() < affectedFactorsVarIndex.size()) {
      // Closest to Root are directly involved keys
      for (auto key : result->markedKeys) {
        constraintGroups[key] = result->markedKeys.size() < updateInvolvedKeys.size() ? 2 : 1;
      }
      // Next are update involved keys
      for (auto key : updateInvolvedKeys) {
        if (!result->markedKeys.count(key)) constraintGroups[key] = 1;
      }
      // After are all convex involved keys and relin keys
    }
    order = Ordering::ColamdConstrained(affectedFactorsVarIndex, constraintGroups);
  } else {
    if (updateParams.constrainedKeys) {
      order = Ordering::ColamdConstrained(affectedFactorsVarIndex, *updateParams.constrainedKeys);
    } else {
      if (theta_.size() > result->observedKeys.size()) {
        // Only if some variables are unconstrained
        FastMap<Key, int> constraintGroups;
        for (Key var : result->observedKeys) constraintGroups[var] = 1;
        order = Ordering::ColamdConstrained(affectedFactorsVarIndex, constraintGroups);
      } else {
        order = Ordering::Colamd(affectedFactorsVarIndex);
      }
    }
  }

  /** (riSAM)
   * On a batch update we still differentiate between factors that are involved because of the update and because they
   * are not convex yet. However, we mark all factors for relin such that we dont use any cached linear factors
   */
  GaussianFactorGraph linearized;
  KeySet relinKeys;
  auto allKeysVec = theta_.keys();
  relinKeys.insert(allKeysVec.begin(), allKeysVec.end());
  for (size_t factor_idx = 0; factor_idx < nonlinearFactors_.size(); factor_idx++) {
    if (RISAM_DEBUG) std::cout << "batch relin factor: " << factor_idx;
    linearized.push_back(relinearizeFactor(factor_idx, updateInvolvedKeys, relinKeys));
    if (RISAM_DEBUG) std::cout << std::endl;
  }

  ISAM2BayesTree::shared_ptr bayesTree =
      ISAM2JunctionTree(GaussianEliminationTree(linearized, affectedFactorsVarIndex, order))
          .eliminate(params_.getEliminationFunction())
          .first;

  roots_.clear();
  roots_.insert(roots_.end(), bayesTree->roots().begin(), bayesTree->roots().end());
  nodes_.clear();
  nodes_.insert(bayesTree->nodes().begin(), bayesTree->nodes().end());

  result->variablesReeliminated = affectedKeysSet->size();
  result->factorsRecalculated = nonlinearFactors_.size();

  // Reeliminated keys for detailed results
  if (params_.enableDetailedResults) {
    for (Key key : theta_.keys()) {
      result->detail->variableStatus[key].isReeliminated = true;
    }
  }
}

/* ************************************************************************* */
void RISAM2::recalculateIncremental(const ISAM2UpdateParams& updateParams, const KeySet& relinKeys,
                                    const FastList<Key>& affectedKeys, const KeySet& updateInvolvedKeys,
                                    KeySet* affectedKeysSet, Cliques* orphans, ISAM2Result* result) {
  const bool debug = ISDEBUG("ISAM2 recalculate");

  // 2. Add the new factors \Factors' into the resulting factor graph
  FastList<Key> affectedAndNewKeys;
  affectedAndNewKeys.insert(affectedAndNewKeys.end(), affectedKeys.begin(), affectedKeys.end());
  affectedAndNewKeys.insert(affectedAndNewKeys.end(), result->observedKeys.begin(), result->observedKeys.end());
  GaussianFactorGraph factors =
      relinearizeAffectedFactors(updateParams, affectedAndNewKeys, relinKeys, updateInvolvedKeys);

  if (debug) {
    factors.print("Relinearized factors: ");
    std::cout << "Affected keys: ";
    for (const Key key : affectedKeys) {
      std::cout << key << " ";
    }
    std::cout << std::endl;
  }

  // Reeliminated keys for detailed results
  if (params_.enableDetailedResults) {
    for (Key key : affectedAndNewKeys) {
      result->detail->variableStatus[key].isReeliminated = true;
    }
  }

  result->variablesReeliminated = affectedAndNewKeys.size();
  result->factorsRecalculated = factors.size();

  // Add the cached intermediate results from the boundary of the orphans...
  GaussianFactorGraph cachedBoundary = UpdateImpl::GetCachedBoundaryFactors(*orphans);
  if (debug) cachedBoundary.print("Boundary factors: ");
  factors.push_back(cachedBoundary);

  // Add the orphaned subtrees
  for (const auto& orphan : *orphans) factors += boost::make_shared<BayesTreeOrphanWrapper<RISAM2::Clique>>(orphan);

  // 3. Re-order and eliminate the factor graph into a Bayes net (Algorithm
  // [alg:eliminate]), and re-assemble into a new Bayes tree (Algorithm
  // [alg:BayesTree])

  // create a partial reordering for the new and contaminated factors
  // result->markedKeys are passed in: those variables will be forced to the
  // end in the ordering
  affectedKeysSet->insert(result->markedKeys.begin(), result->markedKeys.end());
  affectedKeysSet->insert(affectedKeys.begin(), affectedKeys.end());

  VariableIndex affectedFactorsVarIndex(factors);
  FastMap<Key, int> constraintGroups;
  if (risam_params_.use_custom_ordering_constraint) {
    if (updateInvolvedKeys.size() < affectedFactorsVarIndex.size()) {
      // Closest to Root are directly involved keys
      for (auto key : result->markedKeys) {
        constraintGroups[key] = result->markedKeys.size() < updateInvolvedKeys.size() ? 2 : 1;
      }
      // Next are update involved keys
      for (auto key : updateInvolvedKeys) {
        if (!result->markedKeys.count(key)) constraintGroups[key] = 1;
      }
      // After are all convex involved keys and relin keys
    }
  } else {
    if (updateParams.constrainedKeys) {
      constraintGroups = *updateParams.constrainedKeys;
    } else {
      constraintGroups = FastMap<Key, int>();
      const int group = result->observedKeys.size() < affectedFactorsVarIndex.size() ? 1 : 0;
      for (Key var : result->observedKeys) constraintGroups.insert(std::make_pair(var, group));
    }
  }
  // Remove unaffected keys from the constraints
  for (FastMap<Key, int>::iterator iter = constraintGroups.begin(); iter != constraintGroups.end();
       /*Incremented in loop ++iter*/) {
    if (result->unusedKeys.exists(iter->first) || !affectedKeysSet->exists(iter->first))
      constraintGroups.erase(iter++);
    else
      ++iter;
  }

  // Generate ordering
  const Ordering ordering = Ordering::ColamdConstrained(affectedFactorsVarIndex, constraintGroups);

  // Do elimination
  GaussianEliminationTree etree(factors, affectedFactorsVarIndex, ordering);
  auto bayesTree = ISAM2JunctionTree(etree).eliminate(params_.getEliminationFunction()).first;

  roots_.insert(roots_.end(), bayesTree->roots().begin(), bayesTree->roots().end());
  nodes_.insert(bayesTree->nodes().begin(), bayesTree->nodes().end());

  // 4. The orphans have already been inserted during elimination
}

// ##############################################################################################
// ##############################################################################################
// ##############################################################################################
// ##############################################################################################

/* ********************************mode 0(riSAM original)***************************************** */
GaussianFactor::shared_ptr RISAM2::relinearizeFactor_m0(const size_t factor_idx, const KeySet& updateInvolvedKeys,
                                                        const KeySet& relinKeys) {
  auto factor = nonlinearFactors_[factor_idx];
  auto robustFactor = boost::dynamic_pointer_cast<GraduatedFactor>(factor);
  if (!factor) return GaussianFactor::shared_ptr();

  /** Set Inclusion **/

  bool known_inlier = known_inliers_.count(factor_idx);
  bool gnc_factor = !known_inlier;
  bool convex = convex_factors_.count(factor_idx);
  bool update_involved = false;
  bool needs_relin = false;
  for (auto& key : factor->keys()) {
    if (updateInvolvedKeys.find(key) != updateInvolvedKeys.end()) update_involved = true;
    if (relinKeys.find(key) != relinKeys.end()) needs_relin = true;
  }

  /** Determine the type of re-linearization we need to preform **/
  // We use cached linear factors (if configured) for (converged or known-inliers) that dont need relinearzation
  bool update_mu = convex || (gnc_factor && update_involved);
  bool relin = needs_relin || update_mu || !params_.cacheLinearizedFactors;

  /** Update Mu if needed **/
  if (update_mu) {
    if (!robustFactor) {
      throw std::runtime_error("RISAM2 provided a non-graduated factor that is not marked as a known inlier!");
    }

    if (update_involved && new_factors_include_gnc_) {
      factors_to_check_status_.insert(factor_idx);
      mu_[factor_idx] = mu_inits_[factor_idx];
      convex_factors_.insert(factor_idx);
      if (RISAM_DEBUG) std::cout << " | Reset Mu | ";
    } else {
      if (robustFactor->kernel()->isMuConverged(mu_[factor_idx])) convex_factors_.erase(factor_idx);
      mu_[factor_idx] = robustFactor->kernel()->updateMu(mu_[factor_idx]);
      if (RISAM_DEBUG) std::cout << " | Update Mu (" << mu_[factor_idx] << ") | ";
    }
    last_update_info_.mu[factor_idx] = mu_[factor_idx];
  }

  // 2. Linearize the Factor
  GaussianFactor::shared_ptr linearFactor;
  if (relin) {
    if (gnc_factor) {
      linearFactor = robustFactor->linearizeRobust(theta_, mu_[factor_idx]);
      last_update_info_.relinearization_types[factor_idx] = RelinType::ROBUST;
      if (RISAM_DEBUG) std::cout << " | Lin: robust | ";
    } else {
      linearFactor = factor->linearize(theta_);
      last_update_info_.relinearization_types[factor_idx] = RelinType::QUADRATIC;
      if (RISAM_DEBUG) std::cout << " | Lin: Quad | ";
    }
  } else {
    if (RISAM_DEBUG) std::cout << " | Lin: cached | ";
    last_update_info_.relinearization_types[factor_idx] = RelinType::CACHED;
    linearFactor = linearFactors_[factor_idx];
  }

  // Cache regardless of how we linearize
  if (params_.cacheLinearizedFactors) linearFactors_[factor_idx] = linearFactor;
  return linearFactor;
}

/* ********************************B-spline code***************************************** */

double RISAM2::N(int i, int k, double t, const std::vector<double>& knots) {
  /** @brief Obtain basis of B-spline
   * @param i: index of control point
   * @param k: degree
   * @param t: time point (u on paper)
   * @param knots: knots
   */
  int n = knots.size() - 1;
  // Initialize a 2D table to store basis function values
  std::vector<std::vector<double>> basis(n, std::vector<double>(k + 1, 0.0));

  // Compute basis functions for k = 0
  for (int j = 0; j < n; ++j) {
    if (knots[j] <= t && t < knots[j + 1]) {
      basis[j][0] = 1.0;
    } else if (t == knots.back() && j == knots.size() - k - 2) {
      basis[j][0] = 1.0;
    }
  }

  // Compute basis functions for k > 0 using dynamic programming
  for (int d = 1; d <= k; ++d) {
    for (int j = 0; j < n - d; ++j) {
      double term1 = 0.0;
      if (knots[j + d] - knots[j] != 0) {
        term1 = (t - knots[j]) / (knots[j + d] - knots[j]) * basis[j][d - 1];
      }

      double term2 = 0.0;
      if (knots[j + d + 1] - knots[j + 1] != 0) {
        term2 = (knots[j + d + 1] - t) / (knots[j + d + 1] - knots[j + 1]) * basis[j + 1][d - 1];
      }
      basis[j][d] = term1 + term2;
    }
  }

  return basis[i][k];
}

std::vector<Point> RISAM2::bspline_curve(const std::vector<double>& t_values, int degree,
                                         const std::vector<double>& knots, const std::vector<Point>& control_points) {
  /** @brief obtain curve points on B-spline
   * @param t_values: vector of time points 
   * @param degree: degree
   * @param knots: knots
   * @param control_points: control points
   */
  std::vector<Point> curve_points;

  for (double t : t_values) {
    Point curve_point = {0.0, 0.0};

    for (int i = 0; i < control_points.size(); ++i) {
      curve_point.first += N(i, degree, t, knots) * control_points[i].first;
      curve_point.second += N(i, degree, t, knots) * control_points[i].second;
    }
    curve_points.push_back(curve_point);
  }

  return curve_points;
}

Point RISAM2::interpolate(const Point& p1, const Point& p2, double t) {
  return {p1.first + t * (p2.first - p1.first), p1.second + t * (p2.second - p1.second)};
}

std::vector<Point> RISAM2::getInterpolatedControlPoints(double t) {
  std::vector<Point> interpolated;

  for (size_t i = 0; i < control_points_upper.size(); i++) {
    interpolated.push_back(interpolate(control_points_upper[i], control_points_lower[i], t));
  }

  return interpolated;
}

/* ********************************mode 1***************************************** */
GaussianFactor::shared_ptr RISAM2::relinearizeFactor_m1(const size_t factor_idx, const KeySet& updateInvolvedKeys,
                                                        const KeySet& relinKeys) {
  auto factor = nonlinearFactors_[factor_idx];
  auto robustFactor = boost::dynamic_pointer_cast<GraduatedFactor>(factor);
  if (!factor) return GaussianFactor::shared_ptr();

  /** Set Inclusion **/

  bool known_inlier = known_inliers_.count(factor_idx);
  bool gnc_factor = !known_inlier;
  bool convex = convex_factors_.count(factor_idx);
  bool update_involved = false;
  bool needs_relin = false;
  for (auto& key : factor->keys()) {
    if (updateInvolvedKeys.find(key) != updateInvolvedKeys.end()) update_involved = true;
    if (relinKeys.find(key) != relinKeys.end()) needs_relin = true;
  }

  /** Determine the type of re-linearization we need to preform **/
  // We use cached linear factors (if configured) for (converged or known-inliers) that dont need relinearzation
  bool update_mu = convex || (gnc_factor && update_involved);
  bool relin = needs_relin || update_mu || !params_.cacheLinearizedFactors;

  /** Update Mu if needed **/
  if (update_mu) {
    if (!robustFactor) {
      throw std::runtime_error("RISAM2 provided a non-graduated factor that is not marked as a known inlier!");
    }

    if (update_involved && new_factors_include_gnc_) {
      factors_to_check_status_.insert(factor_idx);

      auto residual = robustFactor->residual(theta_);
      auto mahdist = residual;

      auto dist = boost::math::chi_squared_distribution<double>(nonlinearFactors_[factor_idx]->dim());

      double lower_bound = boost::math::quantile(dist, 0.25);
      double upper_bound = boost::math::quantile(dist, 0.90);

      double mu_new = -1;

      num_visited[factor_idx] = init_num_visited[factor_idx];

      if (mahdist < lower_bound) {  // inlier
        std::vector<double> t = {std::min(1.0 / (lower_step_size - 1) * num_visited[factor_idx]++, 1.0)};

        std::vector<Point> result = bspline_curve(t, degree, knots, control_points_middle);

        mu_new = result[0].second;
      } else if (mahdist > upper_bound) {  // outlier
        std::vector<double> t = {std::min(1.0 / (upper_step_size - 1) * num_visited[factor_idx]++, 1.0)};

        std::vector<Point> result = bspline_curve(t, degree, knots, control_points_middle);

        mu_new = result[0].second;
      } else {  // the others
        std::vector<double> t = {std::min(1.0 / (other_step_size - 1) * num_visited[factor_idx]++, 1.0)};

        std::vector<Point> result = bspline_curve(t, degree, knots, control_points_middle);

        mu_new = result[0].second;
      }

      // just in case
      if (mu_new > 1) {
        mu_new = 1;
      }

      mu_[factor_idx] = mu_new;
      convex_factors_.insert(factor_idx);
      if (RISAM_DEBUG) std::cout << " | Reset Mu | ";
    } else {
      if (robustFactor->kernel()->isMuConverged(mu_[factor_idx])) {
        convex_factors_.erase(factor_idx);
      }

      auto residual = robustFactor->residual(theta_);
      auto mahdist = residual;
      auto dist = boost::math::chi_squared_distribution<double>(nonlinearFactors_[factor_idx]->dim());

      double lower_bound = boost::math::quantile(dist, 0.25);
      double upper_bound = boost::math::quantile(dist, 0.90);

      double mu_new = -1;

      if (mahdist < lower_bound) {  // inlier
        std::vector<double> t = {std::min(1.0 / (lower_step_size - 1) * num_visited[factor_idx]++, 1.0)};

        std::vector<Point> result = bspline_curve(t, degree, knots, control_points_middle);

        mu_new = result[0].second;
      } else if (mahdist > upper_bound) {  // outlier
        std::vector<double> t = {std::min(1.0 / (upper_step_size - 1) * num_visited[factor_idx]++, 1.0)};

        std::vector<Point> result = bspline_curve(t, degree, knots, control_points_middle);

        mu_new = result[0].second;
      } else {  // the others
        std::vector<double> t = {std::min(1.0 / (other_step_size - 1) * num_visited[factor_idx]++, 1.0)};

        std::vector<Point> result = bspline_curve(t, degree, knots, control_points_middle);

        mu_new = result[0].second;
      }

      // just in case
      if (mu_new > 1) {
        mu_new = 1;
      }

      mu_[factor_idx] = mu_new;

      if (RISAM_DEBUG) {
        std::cout << " | Update Mu (" << mu_[factor_idx] << ") | ";
      }
    }
    last_update_info_.mu[factor_idx] = mu_[factor_idx];
  }
  // 2. Linearize the Factor
  GaussianFactor::shared_ptr linearFactor;
  if (relin) {
    if (gnc_factor) {
      linearFactor = robustFactor->linearizeRobust(theta_, mu_[factor_idx]);
      last_update_info_.relinearization_types[factor_idx] = RelinType::ROBUST;
      if (RISAM_DEBUG) std::cout << " | Lin: robust | ";
    } else {
      linearFactor = factor->linearize(theta_);
      last_update_info_.relinearization_types[factor_idx] = RelinType::QUADRATIC;
      if (RISAM_DEBUG) std::cout << " | Lin: Quad | ";
    }
  } else {
    if (RISAM_DEBUG) std::cout << " | Lin: cached | ";
    last_update_info_.relinearization_types[factor_idx] = RelinType::CACHED;
    linearFactor = linearFactors_[factor_idx];
  }

  // Cache regardless of how we linearize
  if (params_.cacheLinearizedFactors) linearFactors_[factor_idx] = linearFactor;
  return linearFactor;
}

/* ********************************mode 2***************************************** */
GaussianFactor::shared_ptr RISAM2::relinearizeFactor_m2(const size_t factor_idx, const KeySet& updateInvolvedKeys,
                                                        const KeySet& relinKeys) {
  auto factor = nonlinearFactors_[factor_idx];
  auto robustFactor = boost::dynamic_pointer_cast<GraduatedFactor>(factor);
  if (!factor) return GaussianFactor::shared_ptr();

  /** Set Inclusion **/

  bool known_inlier = known_inliers_.count(factor_idx);
  bool gnc_factor = !known_inlier;
  bool convex = convex_factors_.count(factor_idx);
  bool update_involved = false;
  bool needs_relin = false;
  for (auto& key : factor->keys()) {
    if (updateInvolvedKeys.find(key) != updateInvolvedKeys.end()) update_involved = true;
    if (relinKeys.find(key) != relinKeys.end()) needs_relin = true;
  }

  /** Determine the type of re-linearization we need to preform **/
  // We use cached linear factors (if configured) for (converged or known-inliers) that dont need relinearzation
  bool update_mu = convex || (gnc_factor && update_involved);
  bool relin = needs_relin || update_mu || !params_.cacheLinearizedFactors;

  /** Update Mu if needed **/
  if (update_mu) {
    if (!robustFactor) {
      throw std::runtime_error("RISAM2 provided a non-graduated factor that is not marked as a known inlier!");
    }

    if (update_involved && new_factors_include_gnc_) {
      factors_to_check_status_.insert(factor_idx);

      auto residual = robustFactor->residual(theta_);
      auto mahdist = residual;

      auto dist = boost::math::chi_squared_distribution<double>(nonlinearFactors_[factor_idx]->dim());

      double lower_bound = boost::math::quantile(dist, 0.25);
      double upper_bound = boost::math::quantile(dist, 0.90);

      double mu_new = -1;

      num_visited[factor_idx] = init_num_visited[factor_idx];

      if (mahdist < lower_bound) { // inlier
        control_points = control_points_lower;

        std::vector<double> t = {std::min(1.0 / (lower_step_size - 1) * num_visited[factor_idx]++, 1.0)};

        std::vector<Point> result = bspline_curve(t, degree, knots, control_points);

        mu_new = result[0].second;
      } else if (mahdist > upper_bound) { // outlier
        control_points = control_points_upper;

        std::vector<double> t = {std::min(1.0 / (upper_step_size - 1) * num_visited[factor_idx]++, 1.0)};

        std::vector<Point> result = bspline_curve(t, degree, knots, control_points);

        mu_new = result[0].second;
      } else { // the others
        double prob = boost::math::cdf(dist, mahdist);
        prob = (prob - 0.25) / (0.9 - 0.25);
        control_points = getInterpolatedControlPoints(prob);

        std::vector<double> t = {std::min(1.0 / (other_step_size - 1) * num_visited[factor_idx]++, 1.0)};

        std::vector<Point> result = bspline_curve(t, degree, knots, control_points);

        mu_new = result[0].second;
      }

      // just in case
      if (mu_new > 1) {
        mu_new = 1;
      }

      mu_[factor_idx] = mu_new;
      convex_factors_.insert(factor_idx);
      if (RISAM_DEBUG) std::cout << " | Reset Mu | ";
    } else {
      if (robustFactor->kernel()->isMuConverged(mu_[factor_idx])) {
        convex_factors_.erase(factor_idx);
      }

      auto residual = robustFactor->residual(theta_);
      auto mahdist = residual;
      auto dist = boost::math::chi_squared_distribution<double>(nonlinearFactors_[factor_idx]->dim());

      double lower_bound = boost::math::quantile(dist, 0.25);
      double upper_bound = boost::math::quantile(dist, 0.90);

      double mu_new = -1;

      if (mahdist < lower_bound) { // inlier
        control_points = control_points_lower;

        std::vector<double> t = {std::min(1.0 / (lower_step_size - 1) * num_visited[factor_idx]++, 1.0)};

        std::vector<Point> result = bspline_curve(t, degree, knots, control_points);

        mu_new = result[0].second;
      } else if (mahdist > upper_bound) { // outlier
        control_points = control_points_upper;

        std::vector<double> t = {std::min(1.0 / (upper_step_size - 1) * num_visited[factor_idx]++, 1.0)};

        std::vector<Point> result = bspline_curve(t, degree, knots, control_points);

        mu_new = result[0].second;
      } else { // the others
        double prob = boost::math::cdf(dist, mahdist);
        prob = (prob - 0.25) / (0.9 - 0.25);
        control_points = getInterpolatedControlPoints(prob);

        std::vector<double> t = {std::min(1.0 / (other_step_size - 1) * num_visited[factor_idx]++, 1.0)};

        std::vector<Point> result = bspline_curve(t, degree, knots, control_points);

        mu_new = result[0].second;
      }

      // just in case
      if (mu_new > 1) {
        mu_new = 1;
      }

      mu_[factor_idx] = mu_new;

      if (RISAM_DEBUG) {
        std::cout << " | Update Mu (" << mu_[factor_idx] << ") | ";
      }
    }
    last_update_info_.mu[factor_idx] = mu_[factor_idx];
  }

  // 2. Linearize the Factor
  GaussianFactor::shared_ptr linearFactor;
  if (relin) {
    if (gnc_factor) {
      linearFactor = robustFactor->linearizeRobust(theta_, mu_[factor_idx]);
      last_update_info_.relinearization_types[factor_idx] = RelinType::ROBUST;
      if (RISAM_DEBUG) std::cout << " | Lin: robust | ";
    } else {
      linearFactor = factor->linearize(theta_);
      last_update_info_.relinearization_types[factor_idx] = RelinType::QUADRATIC;
      if (RISAM_DEBUG) std::cout << " | Lin: Quad | ";
    }
  } else {
    if (RISAM_DEBUG) std::cout << " | Lin: cached | ";
    last_update_info_.relinearization_types[factor_idx] = RelinType::CACHED;
    linearFactor = linearFactors_[factor_idx];
  }

  // Cache regardless of how we linearize
  if (params_.cacheLinearizedFactors) linearFactors_[factor_idx] = linearFactor;
  return linearFactor;
}


/* ********************************find optimal mu***************************************** */
double RISAM2::partial_derivative(double mu, double r, int c) {
    double term1 = (4 * std::pow(c, 2) * std::pow(mu, 2) * std::pow(r, 4 * mu)) / std::pow(std::pow(c, 2) + std::pow(r, 2 * mu), 3);
    double term2 = (2 * std::pow(c, 2) * std::pow(mu, 2) * std::pow(r, 2 * mu)) / std::pow(std::pow(c, 2) + std::pow(r, 2 * mu), 2);
    double term3 = (3 * std::pow(c, 2) * mu * std::pow(r, 2 * mu)) / std::pow(std::pow(c, 2) + std::pow(r, 2 * mu), 2);
    double term4 = std::pow(c, 2) / (std::pow(c, 2) + std::pow(r, 2 * mu));
    return term1 - term2 - term3 + term4;
}

double RISAM2::bisection(double a, double b, double tol, double r, int dim) {
    double fa = partial_derivative(a, r, dim);
    double fb = partial_derivative(b, r, dim);

    if (fa * fb >= 0) {
        std::cout << "Bisection method fails." << std::endl;
        return -1;
    }

    double c = a;
    while ((b - a) / 2.0 > tol) {
        c = (a + b) / 2.0;
        double fc = partial_derivative(c, r, dim);

        if (fc == 0) {
            break;
        } else if (fc * fa < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    return c;
}


/* ********************************mode 3***************************************** */
GaussianFactor::shared_ptr RISAM2::relinearizeFactor_m3(const size_t factor_idx, const KeySet& updateInvolvedKeys,
                                                        const KeySet& relinKeys) {
  auto factor = nonlinearFactors_[factor_idx];
  auto robustFactor = boost::dynamic_pointer_cast<GraduatedFactor>(factor);
  if (!factor) return GaussianFactor::shared_ptr();

  /** Set Inclusion **/

  bool known_inlier = known_inliers_.count(factor_idx);
  bool gnc_factor = !known_inlier;
  bool convex = convex_factors_.count(factor_idx);
  bool update_involved = false;
  bool needs_relin = false;
  for (auto& key : factor->keys()) {
    if (updateInvolvedKeys.find(key) != updateInvolvedKeys.end()) update_involved = true;
    if (relinKeys.find(key) != relinKeys.end()) needs_relin = true;
  }

  /** Determine the type of re-linearization we need to preform **/
  // We use cached linear factors (if configured) for (converged or known-inliers) that dont need relinearzation
  bool update_mu = convex || (gnc_factor && update_involved);
  bool relin = needs_relin || update_mu || !params_.cacheLinearizedFactors;

  double lower = 0.0;
  double upper = 1.0;
  double tol = 1e-9;

  /** Update Mu if needed **/
  if (update_mu) {
    if (!robustFactor) {
      throw std::runtime_error("RISAM2 provided a non-graduated factor that is not marked as a known inlier!");
    }

    if (update_involved && new_factors_include_gnc_) {
      factors_to_check_status_.insert(factor_idx);

      auto residual = robustFactor->residual(theta_);
      auto mahdist = residual;

      auto dist = boost::math::chi_squared_distribution<double>(nonlinearFactors_[factor_idx]->dim());

      double lower_bound = boost::math::quantile(dist, 0.25);
      double upper_bound = boost::math::quantile(dist, 0.90);

      double mu_new = -1;

      num_visited[factor_idx] = init_num_visited[factor_idx];

      if (mahdist > upper_bound || mahdist < sqrt(3)) { // outlier or inlier, only convex optimization.
        if (num_visited[factor_idx] == 0)
          mu_new = 0;
        else
          mu_new = 1;
      } else { // the other, convex(0), change convex to non-convex(mu), non-convex optimization(1)
        if (num_visited[factor_idx] == 0)
          mu_new = 0;
        else if (num_visited[factor_idx] == 1)
        { // no case
          mu_new = bisection(lower, upper, tol, mahdist, nonlinearFactors_[factor_idx]->dim());
          std::cout << mu_new << "********************************** PROBLEM !! **************************************" << std::endl;
        }
        else
          mu_new = 1;
      }

      num_visited[factor_idx]++;

      // just in case
      if (mu_new > 1) {
        std::cout << "just in case" << std::endl;
        mu_new = 1;
      }

      mu_[factor_idx] = mu_new;
      convex_factors_.insert(factor_idx);
      if (RISAM_DEBUG) std::cout << " | Reset Mu | ";
    } else {
      if (robustFactor->kernel()->isMuConverged(mu_[factor_idx])) {
        convex_factors_.erase(factor_idx);
      }

      auto residual = robustFactor->residual(theta_);
      auto mahdist = residual;
      auto dist = boost::math::chi_squared_distribution<double>(nonlinearFactors_[factor_idx]->dim());

      double lower_bound = boost::math::quantile(dist, 0.25);
      double upper_bound = boost::math::quantile(dist, 0.90);

      double mu_new = -1;

      if (mahdist > upper_bound || mahdist < sqrt(3)) { // outlier or inlier, only convex optimization.
        if (num_visited[factor_idx] == 0)
          mu_new = 0;
        else
          mu_new = 1;
      } else { // the other, convex(0), change convex to non-convex(mu), non-convex optimization(1)
        if (num_visited[factor_idx] == 0)
          mu_new = 0;
        else if (num_visited[factor_idx] == 1)
        {
          mu_new = bisection(lower, upper, tol, mahdist, nonlinearFactors_[factor_idx]->dim());
          if (mu_new > 1 || mu_new < 0)
            std::cout << mu_new << "********************************** PROBLEM !! **************************************" << std::endl;
        }
        else
          mu_new = 1;
      }

      num_visited[factor_idx]++;

      // just in case
      if (mu_new > 1) {
        std::cout << "just in case" << std::endl;
        mu_new = 1;
      }

      mu_[factor_idx] = mu_new;

      if (RISAM_DEBUG) {
        std::cout << " | Update Mu (" << mu_[factor_idx] << ") | ";
      }
    }
    last_update_info_.mu[factor_idx] = mu_[factor_idx];
  }

  // 2. Linearize the Factor
  GaussianFactor::shared_ptr linearFactor;
  if (relin) {
    if (gnc_factor) {
      linearFactor = robustFactor->linearizeRobust(theta_, mu_[factor_idx]);
      last_update_info_.relinearization_types[factor_idx] = RelinType::ROBUST;
      if (RISAM_DEBUG) std::cout << " | Lin: robust | ";
    } else {
      linearFactor = factor->linearize(theta_);
      last_update_info_.relinearization_types[factor_idx] = RelinType::QUADRATIC;
      if (RISAM_DEBUG) std::cout << " | Lin: Quad | ";
    }
  } else {
    if (RISAM_DEBUG) std::cout << " | Lin: cached | ";
    last_update_info_.relinearization_types[factor_idx] = RelinType::CACHED;
    linearFactor = linearFactors_[factor_idx];
  }

  // Cache regardless of how we linearize
  if (params_.cacheLinearizedFactors) linearFactors_[factor_idx] = linearFactor;
  return linearFactor;
}

// ##############################################################################################
// ##############################################################################################
// ##############################################################################################
// ##############################################################################################

/* ************************************************************************* */
GaussianFactorGraph RISAM2::relinearizeAffectedFactors(const ISAM2UpdateParams& updateParams,
                                                       const FastList<Key>& affectedKeys, const KeySet& relinKeys,
                                                       const KeySet& updateInvolvedKeys) {
  FactorIndexSet candidates = UpdateImpl::GetAffectedFactors(affectedKeys, variableIndex_);

  // for fast lookup below
  KeySet affectedKeysSet;
  affectedKeysSet.insert(affectedKeys.begin(), affectedKeys.end());

  GaussianFactorGraph linearized;
  for (const FactorIndex idx : candidates) {
    bool inside = true;
    for (Key key : nonlinearFactors_[idx]->keys()) {
      if (affectedKeysSet.find(key) == affectedKeysSet.end()) {
        inside = false;
        break;
      }
    }
    if (inside) {
      if (RISAM_DEBUG) std::cout << "Inc relin factor: " << idx;
      linearized.push_back(relinearizeFactor(idx, updateInvolvedKeys, relinKeys));
      if (RISAM_DEBUG) std::cout << std::endl;
    }
  }

  return linearized;
}

// ##############################################################################################
// ##############################################################################################
// ##############################################################################################
// ##############################################################################################

/* ********************************mode 0***************************************** */
void RISAM2::pushBackFactors_m0(const NonlinearFactorGraph& newFactors, bool known_inliers,
                                NonlinearFactorGraph* nonlinearFactors, GaussianFactorGraph* linearFactors,
                                VariableIndex* variableIndex, FactorIndices* newFactorsIndices,
                                KeySet* keysWithRemovedFactors, ISAM2UpdateParams updateParams) {
  // Perform the first part of the bookkeeping updates for adding new factors.
  // Adds them to the complete list of nonlinear factors, and populates the
  // list of new factor indices, both optionally finding and reusing empty
  // factor slots.
  *newFactorsIndices = nonlinearFactors->add_factors(newFactors, params_.findUnusedFactorSlots);

  for (auto idx : *newFactorsIndices) {
    double mu_init = 0.0;
    auto robustFactor = boost::dynamic_pointer_cast<risam::GraduatedFactor>(nonlinearFactors_[idx]);
    if (robustFactor) {
      mu_init = robustFactor->kernel()->muInit();
    }

    if (idx >= mu_.size()) {
      mu_inits_.push_back(mu_init);
      mu_.push_back(mu_init);
    } else {
      mu_inits_[idx] = 0.0;
      mu_[idx] = 0.0;
    }
    if (known_inliers) known_inliers_.insert(idx);
  }

  // Remove the removed factors
  NonlinearFactorGraph removedFactors;
  removedFactors.reserve(updateParams.removeFactorIndices.size());
  for (const auto index : updateParams.removeFactorIndices) {
    removedFactors.push_back(nonlinearFactors->at(index));
    nonlinearFactors->remove(index);
    if (params_.cacheLinearizedFactors) linearFactors->remove(index);
  }

  // Remove removed factors from the variable index so we do not attempt to
  // relinearize them
  variableIndex->remove(updateParams.removeFactorIndices.begin(), updateParams.removeFactorIndices.end(),
                        removedFactors);
  *keysWithRemovedFactors = removedFactors.keys();
}

/* ********************************mode 1***************************************** */
void RISAM2::pushBackFactors_m1(const NonlinearFactorGraph& newFactors, bool known_inliers,
                                NonlinearFactorGraph* nonlinearFactors, GaussianFactorGraph* linearFactors,
                                VariableIndex* variableIndex, FactorIndices* newFactorsIndices,
                                KeySet* keysWithRemovedFactors, ISAM2UpdateParams updateParams) {
  // Perform the first part of the bookkeeping updates for adding new factors.
  // Adds them to the complete list of nonlinear factors, and populates the
  // list of new factor indices, both optionally finding and reusing empty
  // factor slots.
  *newFactorsIndices = nonlinearFactors->add_factors(newFactors, params_.findUnusedFactorSlots);

  for (auto idx : *newFactorsIndices) {
    double mu_init = 0.0;
    auto robustFactor = boost::dynamic_pointer_cast<risam::GraduatedFactor>(nonlinearFactors_[idx]);
    if (robustFactor) {
      mu_init = robustFactor->kernel()->muInit();
    }

    if (idx >= mu_.size()) {
      init_num_visited[idx] = 0;
      mu_inits_.push_back(mu_init);
      mu_.push_back(mu_init);
    } else {
      init_num_visited[idx] = 0;
      mu_inits_[idx] = 0.0;
      mu_[idx] = 0.0;
    }
    if (known_inliers) known_inliers_.insert(idx);
  }

  // Remove the removed factors
  NonlinearFactorGraph removedFactors;
  removedFactors.reserve(updateParams.removeFactorIndices.size());
  for (const auto index : updateParams.removeFactorIndices) {
    removedFactors.push_back(nonlinearFactors->at(index));
    nonlinearFactors->remove(index);
    if (params_.cacheLinearizedFactors) linearFactors->remove(index);
  }

  // Remove removed factors from the variable index so we do not attempt to
  // relinearize them
  variableIndex->remove(updateParams.removeFactorIndices.begin(), updateParams.removeFactorIndices.end(),
                        removedFactors);
  *keysWithRemovedFactors = removedFactors.keys();
}

/* ********************************mode 2***************************************** */
void RISAM2::pushBackFactors_m2(const NonlinearFactorGraph& newFactors, bool known_inliers,
                                NonlinearFactorGraph* nonlinearFactors, GaussianFactorGraph* linearFactors,
                                VariableIndex* variableIndex, FactorIndices* newFactorsIndices,
                                KeySet* keysWithRemovedFactors, ISAM2UpdateParams updateParams) {
  // Perform the first part of the bookkeeping updates for adding new factors.
  // Adds them to the complete list of nonlinear factors, and populates the
  // list of new factor indices, both optionally finding and reusing empty
  // factor slots.
  *newFactorsIndices = nonlinearFactors->add_factors(newFactors, params_.findUnusedFactorSlots);

  for (auto idx : *newFactorsIndices) {
    double mu_init = 0.0;
    auto robustFactor = boost::dynamic_pointer_cast<risam::GraduatedFactor>(nonlinearFactors_[idx]);
    if (robustFactor) {
      mu_init = robustFactor->kernel()->muInit();
    }

    if (idx >= mu_.size()) {
      init_num_visited[idx] = 0;
      mu_inits_.push_back(mu_init);
      mu_.push_back(mu_init);
    } else {
      init_num_visited[idx] = 0;
      mu_inits_[idx] = 0.0;
      mu_[idx] = 0.0;
    }
    if (known_inliers) known_inliers_.insert(idx);
  }

  // Remove the removed factors
  NonlinearFactorGraph removedFactors;
  removedFactors.reserve(updateParams.removeFactorIndices.size());
  for (const auto index : updateParams.removeFactorIndices) {
    removedFactors.push_back(nonlinearFactors->at(index));
    nonlinearFactors->remove(index);
    if (params_.cacheLinearizedFactors) linearFactors->remove(index);
  }

  // Remove removed factors from the variable index so we do not attempt to
  // relinearize them
  variableIndex->remove(updateParams.removeFactorIndices.begin(), updateParams.removeFactorIndices.end(),
                        removedFactors);
  *keysWithRemovedFactors = removedFactors.keys();
}

/* ********************************mode 3***************************************** */
void RISAM2::pushBackFactors_m3(const NonlinearFactorGraph& newFactors, bool known_inliers,
                                NonlinearFactorGraph* nonlinearFactors, GaussianFactorGraph* linearFactors,
                                VariableIndex* variableIndex, FactorIndices* newFactorsIndices,
                                KeySet* keysWithRemovedFactors, ISAM2UpdateParams updateParams) {
  // Perform the first part of the bookkeeping updates for adding new factors.
  // Adds them to the complete list of nonlinear factors, and populates the
  // list of new factor indices, both optionally finding and reusing empty
  // factor slots.
  *newFactorsIndices = nonlinearFactors->add_factors(newFactors, params_.findUnusedFactorSlots);

  for (auto idx : *newFactorsIndices) {
    double mu_init = 0.0;
    auto robustFactor = boost::dynamic_pointer_cast<risam::GraduatedFactor>(nonlinearFactors_[idx]);
    if (robustFactor) {
      mu_init = robustFactor->kernel()->muInit();
    }

    if (idx >= mu_.size()) {
      init_num_visited[idx] = 0;
      mu_inits_.push_back(mu_init);
      mu_.push_back(mu_init);
    } else {
      init_num_visited[idx] = 0;
      mu_inits_[idx] = 0.0;
      mu_[idx] = 0.0;
    }
    if (known_inliers) known_inliers_.insert(idx);
  }

  // Remove the removed factors
  NonlinearFactorGraph removedFactors;
  removedFactors.reserve(updateParams.removeFactorIndices.size());
  for (const auto index : updateParams.removeFactorIndices) {
    removedFactors.push_back(nonlinearFactors->at(index));
    nonlinearFactors->remove(index);
    if (params_.cacheLinearizedFactors) linearFactors->remove(index);
  }

  // Remove removed factors from the variable index so we do not attempt to
  // relinearize them
  variableIndex->remove(updateParams.removeFactorIndices.begin(), updateParams.removeFactorIndices.end(),
                        removedFactors);
  *keysWithRemovedFactors = removedFactors.keys();
}

// ##############################################################################################
// ##############################################################################################
// ##############################################################################################
// ##############################################################################################

/* ************************************************************************* */
double RISAM2::robustError(boost::optional<VectorValues> delta) const {
  Values theta = theta_;
  if (delta) {
    theta = theta.retract(*delta);
  }
  return robustError(theta);
}

/* ************************************************************************* */
double RISAM2::robustError(Values vals) const {  // TODO (dan) Optimize
  double error = 0;
  for (size_t i = 0; i < nonlinearFactors_.size(); i++) {
    auto factor = nonlinearFactors_[i];
    auto robustFactor = boost::dynamic_pointer_cast<risam::GraduatedFactor>(nonlinearFactors_[i]);
    if (robustFactor) {
      error += robustFactor->robustResidual(vals, mu_[i]);
    } else if (factor) {
      error += sqrt(2.0 * factor->error(vals));
    }
  }
  return error;
}

/**
 * ########  ########    ######## ########     ###    ##     ## ######## ########   ######     ###    ##
 * ##     ##    ##          ##    ##     ##   ## ##   ##     ## ##       ##     ## ##    ##   ## ##   ##
 * ##     ##    ##          ##    ##     ##  ##   ##  ##     ## ##       ##     ## ##        ##   ##  ##
 * ########     ##          ##    ########  ##     ## ##     ## ######   ########   ######  ##     ## ##
 * ##     ##    ##          ##    ##   ##   #########  ##   ##  ##       ##   ##         ## ######### ##
 * ##     ##    ##          ##    ##    ##  ##     ##   ## ##   ##       ##    ##  ##    ## ##     ## ##
 * ########     ##          ##    ##     ## ##     ##    ###    ######## ##     ##  ######  ##     ## ########
 */

KeySet RISAM2::traverseTop(const KeyVector& keys) {
  KeySet traversedKeys;
  // process each key of the new factor
  for (const Key& j : keys) {
    typename Nodes::const_iterator node = nodes_.find(j);
    if (node != nodes_.end()) {
      // remove path from clique to root
      this->traversePath(traversedKeys, node->second);
    }
  }
  return traversedKeys;
}

void RISAM2::traversePath(KeySet& traversedKeys, ISAM2::sharedClique clique) {
  // base case is nullptr, if so we do nothing and return empties above
  if (clique) {
    // traverse me
    this->traverseClique(traversedKeys, clique);
    // traverse path above me
    this->traversePath(traversedKeys, typename ISAM2::Clique::shared_ptr(clique->parent_.lock()));
  }
}

void RISAM2::traverseClique(KeySet& traversedKeys, ISAM2::sharedClique clique) {
  traversedKeys.insert(clique->conditional()->frontals().begin(), clique->conditional()->frontals().end());
}

/**
 *  ######   #######  ##       ##     ## ########
 * ##    ## ##     ## ##       ##     ## ##
 * ##       ##     ## ##       ##     ## ##
 *  ######  ##     ## ##       ##     ## ######
 *       ## ##     ## ##        ##   ##  ##
 * ##    ## ##     ## ##         ## ##   ##
 *  ######   #######  ########    ###    ########
 */

/* ************************************************************************* */
// Marked const but actually changes mutable delta
void RISAM2::updateDelta(bool forceFullSolve) const {
  if (!risam_params_
           .optimization_params) {  // if RISAM does not provide optimizer params, use the ISAM2 configured ones
    if (params_.optimizationParams.type() == typeid(ISAM2GaussNewtonParams)) {
      // If using Gauss-Newton, update with wildfireThreshold
      const ISAM2GaussNewtonParams& gaussNewtonParams = boost::get<ISAM2GaussNewtonParams>(params_.optimizationParams);
      const double effectiveWildfireThreshold = forceFullSolve ? 0.0 : gaussNewtonParams.wildfireThreshold;

      DeltaImpl::UpdateGaussNewtonDelta(roots_, deltaReplacedMask_, effectiveWildfireThreshold, &delta_);
      deltaReplacedMask_.clear();

    } else if (params_.optimizationParams.type() == typeid(ISAM2DoglegParams)) {
      // If using Dogleg, do a Dogleg step
      const ISAM2DoglegParams& doglegParams = boost::get<ISAM2DoglegParams>(params_.optimizationParams);
      const double effectiveWildfireThreshold = forceFullSolve ? 0.0 : doglegParams.wildfireThreshold;

      // Do one Dogleg iteration

      // Compute Newton's method step
      DeltaImpl::UpdateGaussNewtonDelta(roots_, deltaReplacedMask_, effectiveWildfireThreshold, &deltaNewton_);

      // Compute steepest descent step
      const VectorValues gradAtZero = this->gradientAtZero();
      DeltaImpl::UpdateRgProd(roots_, deltaReplacedMask_, gradAtZero, &RgProd_);
      const VectorValues dx_u = DeltaImpl::ComputeGradientSearch(gradAtZero, RgProd_);
      deltaReplacedMask_.clear();

      // Compute dogleg point
      struct RobustErrorStruct {
        NonlinearFactorGraph graph_;
        FastVector<double> mu_;
        RobustErrorStruct(NonlinearFactorGraph graph, FastVector<double> mu) {
          graph_ = graph;
          mu_ = mu;
        }
        double error(const Values& theta) const {
          double error = 0;
          for (size_t i = 0; i < graph_.size(); i++) {
            auto factor = graph_[i];
            auto robustFactor = boost::dynamic_pointer_cast<risam::GraduatedFactor>(graph_[i]);
            if (robustFactor) {
              error += robustFactor->robustResidual(theta, mu_[i]);
            } else if (factor) {
              error += sqrt(2.0 * factor->error(theta));
            }
          }
          return error;
        }
      };

      auto func = RobustErrorStruct(nonlinearFactors_, mu_);
      DoglegOptimizerImpl::IterationResult doglegResult(
          DoglegOptimizerImpl::Iterate(*doglegDelta_, doglegParams.adaptationMode, dx_u, deltaNewton_, *this, func,
                                       theta_, func.error(theta_), doglegParams.verbose));

      // Update Delta and linear step

      doglegDelta_ = doglegResult.delta;
      delta_ = doglegResult.dx_d;

    } else {
      throw std::runtime_error("iSAM2: unknown ISAM2Params type");
    }
  } else {
    // If using Dogleg, do a Dogleg step
    const DoglegLineSearchParams& linesearch_params = *risam_params_.optimization_params;
    const double effectiveWildfireThreshold = forceFullSolve ? 0.0 : linesearch_params.wildfire_threshold;

    // Do one DoglegLineSearch iteration#

    // Compute Newton's method step

    DeltaImpl::UpdateGaussNewtonDelta(roots_, deltaReplacedMask_, effectiveWildfireThreshold, &deltaNewton_);

    // Compute steepest descent step
    const VectorValues gradAtZero = this->gradientAtZero();
    DeltaImpl::UpdateRgProd(roots_, deltaReplacedMask_, gradAtZero, &RgProd_);
    const VectorValues dx_u = DeltaImpl::ComputeGradientSearch(gradAtZero, RgProd_);
    deltaReplacedMask_.clear();

    // Compute dogleg point
    DoglegOptimizerImpl::IterationResult doglegResult(
        DoglegLineSearchImpl::Iterate(linesearch_params, dx_u, deltaNewton_, *this, theta_));

    // Update Delta and linear step

    delta_ = doglegResult.dx_d;
  }
}

/* ************************************************************************* */
void RISAM2::iterateToConvergence() {
  // Iterate At least until mu parameters have converged
  while (risam_params_.converge_mu && convex_factors_.size()) {
    update(NonlinearFactorGraph(), Values(), false, ISAM2UpdateParams());
  }
  // Iterate until delta update is sufficiently small
  if (risam_params_.converge_values) {
    for (size_t i = 0; i < risam_params_.value_converge_max_iters; i++) {
      double error = robustError(delta_);
      if (error < risam_params_.value_converge_abs_tol) break;
      update(NonlinearFactorGraph(), Values(), false, ISAM2UpdateParams());
      double new_error = robustError(delta_);
      if (error - new_error < risam_params_.value_converge_abs_tol ||
          (error - new_error) / error < risam_params_.value_converge_rel_tol)
        break;
    }
  }
}

// ##############################################################################################
// ##############################################################################################
// ##############################################################################################
// ##############################################################################################

/* ********************************mode 0***************************************** */
void RISAM2::updateMuInitOnConvergence_m0() {
  for (auto factor_idx : factors_to_check_status_) {
    auto robustFactor = boost::dynamic_pointer_cast<GraduatedFactor>(nonlinearFactors_[factor_idx]);
    auto residual = robustFactor->residual(theta_);
    auto mahdist = residual;  // Mahalanobis Distance
    auto dist = boost::math::chi_squared_distribution<double>(nonlinearFactors_[factor_idx]->dim());
    if (mahdist > boost::math::quantile(dist, risam_params_.outlier_mu_chisq_upper_bound)) {
      mu_inits_[factor_idx] = robustFactor->kernel()->updateMu(mu_inits_[factor_idx]);
    } else if (mahdist < boost::math::quantile(dist, risam_params_.outlier_mu_chisq_lower_bound)) {
      mu_inits_[factor_idx] = robustFactor->kernel()->updateMuInv(mu_inits_[factor_idx]);
    }
  }
  factors_to_check_status_.clear();
}

/* ********************************mode 1***************************************** */
void RISAM2::updateMuInitOnConvergence_m1() {
  for (auto factor_idx : factors_to_check_status_) {
    auto robustFactor = boost::dynamic_pointer_cast<GraduatedFactor>(nonlinearFactors_[factor_idx]);
    auto residual = robustFactor->residual(theta_);
    auto mahdist = residual;  // Mahalanobis Distance
    auto dist = boost::math::chi_squared_distribution<double>(nonlinearFactors_[factor_idx]->dim());
    if (mahdist > boost::math::quantile(dist, risam_params_.outlier_mu_chisq_upper_bound)) { // outlier
      init_num_visited[factor_idx] = 10;
    } else if (mahdist < boost::math::quantile(dist, risam_params_.outlier_mu_chisq_lower_bound)) { // inlier
      if (init_num_visited[factor_idx] >= 1) init_num_visited[factor_idx] = 0;
    }
  }
  factors_to_check_status_.clear();
}

/* ********************************mode 2***************************************** */
void RISAM2::updateMuInitOnConvergence_m2() {
  for (auto factor_idx : factors_to_check_status_) {
    auto robustFactor = boost::dynamic_pointer_cast<GraduatedFactor>(nonlinearFactors_[factor_idx]);
    auto residual = robustFactor->residual(theta_);
    auto mahdist = residual;  // Mahalanobis Distance
    auto dist = boost::math::chi_squared_distribution<double>(nonlinearFactors_[factor_idx]->dim());
    if (mahdist > boost::math::quantile(dist, risam_params_.outlier_mu_chisq_upper_bound)) { // outlier
      init_num_visited[factor_idx] = 10;
    } else if (mahdist < boost::math::quantile(dist, risam_params_.outlier_mu_chisq_lower_bound)) { // inlier
      if (init_num_visited[factor_idx] >= 1) init_num_visited[factor_idx] = 0;
    }
  }
  factors_to_check_status_.clear();
}

/* ********************************mode 3***************************************** */
void RISAM2::updateMuInitOnConvergence_m3() {
  for (auto factor_idx : factors_to_check_status_) {
    auto robustFactor = boost::dynamic_pointer_cast<GraduatedFactor>(nonlinearFactors_[factor_idx]);
    auto residual = robustFactor->residual(theta_);
    auto mahdist = residual;  // Mahalanobis Distance
    auto dist = boost::math::chi_squared_distribution<double>(nonlinearFactors_[factor_idx]->dim());
    if (mahdist > boost::math::quantile(dist, risam_params_.outlier_mu_chisq_upper_bound)) { // outlier
      init_num_visited[factor_idx] = 2;
    } else if (mahdist < boost::math::quantile(dist, risam_params_.outlier_mu_chisq_lower_bound)) { // inlier
      if (init_num_visited[factor_idx] >= 1) init_num_visited[factor_idx] = 0;
    }
  }
  factors_to_check_status_.clear();
}

// ##############################################################################################
// ##############################################################################################
// ##############################################################################################
// ##############################################################################################

/* ************************************************************************* */
VectorValues& RISAM2::getDelta() const {
  // Delta is updated after each update, so it is always up-to-date
  return delta_;
}

/* ************************************************************************* */
Values RISAM2::calculateEstimate() const {
  const VectorValues& delta(getDelta());
  return theta_.retract(delta);
}

/**
 * ##     ## ####  ######  ##     ##    ###    ##       #### ######## ########
 * ##     ##  ##  ##    ## ##     ##   ## ##   ##        ##       ##  ##
 * ##     ##  ##  ##       ##     ##  ##   ##  ##        ##      ##   ##
 * ##     ##  ##   ######  ##     ## ##     ## ##        ##     ##    ######
 *  ##   ##   ##        ## ##     ## ######### ##        ##    ##     ##
 *   ## ##    ##  ##    ## ##     ## ##     ## ##        ##   ##      ##
 *    ###    ####  ######   #######  ##     ## ######## #### ######## ########
 */

/* ************************************************************************* */
void RISAM2::dot(std::ostream& os, const KeyFormatter& keyFormatter) const {
  if (roots_.empty()) throw std::invalid_argument("the root of Bayes tree has not been initialized!");
  os << "digraph G{\n";
  for (const sharedClique& root : roots_) dot(os, root, keyFormatter, 0);
  os << "}";
  std::flush(os);
}

/* ************************************************************************* */
std::string RISAM2::dot(const KeyFormatter& keyFormatter) const {
  std::stringstream ss;
  dot(ss, keyFormatter);
  return ss.str();
}

/* ************************************************************************* */
void RISAM2::saveGraph(const std::string& filename, const KeyFormatter& keyFormatter) const {
  std::ofstream of(filename.c_str());
  dot(of, keyFormatter);
  of.close();
}

/* ************************************************************************* */
void RISAM2::dot(std::ostream& s, sharedClique clique, const KeyFormatter& keyFormatter, int parentnum) const {
  static int num = 0;
  bool first = true;
  std::stringstream out;
  out << num;
  std::string parent = out.str();
  parent += "[label=\"";

  for (Key key : clique->conditional_->frontals()) {
    if (!first) parent += ", ";
    first = false;
    parent += keyFormatter(key);
  }

  if (clique->parent()) {
    parent += " : ";
    s << parentnum << "->" << num << "\n";
  }

  first = true;
  for (Key parentKey : clique->conditional_->parents()) {
    if (!first) parent += ", ";
    first = false;
    parent += keyFormatter(parentKey);
  }

  parent += "\"";
  bool update = false;
  bool convex = false;
  for (Key key : clique->conditional_->frontals()) {
    if (last_update_info_.updateInvolvedKeys.find(key) != last_update_info_.updateInvolvedKeys.end()) update = true;
    if (last_update_info_.affectedKeysConvex.find(key) != last_update_info_.affectedKeysConvex.end()) convex = true;
  }
  for (Key key : clique->conditional_->parents()) {
    if (last_update_info_.updateInvolvedKeys.find(key) != last_update_info_.updateInvolvedKeys.end()) update = true;
    if (last_update_info_.affectedKeysConvex.find(key) != last_update_info_.affectedKeysConvex.end()) convex = true;
  }
  if (update) {
    parent += ", color=firebrick";
  } else if (convex) {
    parent += ", color=dodgerblue2";
  }
  parent += "];\n";

  s << parent;
  parentnum = num;

  for (sharedClique c : clique->children) {
    num++;
    dot(s, c, keyFormatter, parentnum);
  }
}

}  // namespace risam