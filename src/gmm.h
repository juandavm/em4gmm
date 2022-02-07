/* Expectation Maximization for Gaussian Mixture Models.
Copyright (C) 2012-2022 Juan Valor

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */

#ifndef _gmm_h
#define _gmm_h

	typedef struct{
		decimal prior;  /* Prior probability of this component of the mixture.   */
		decimal cgauss; /* Cache storage in order to make faster the classifier. */
		decimal *mean;  /* Means vector of gaussian multivariate distribution.   */
		decimal *dcov;  /* Diagonal covariances, when classify are the inverses. */
		decimal *_mean; /* Counts to estimate the future parameter (used by EM). */
		decimal *_dcov; /* Counts to estimate the future parameter (used by EM). */
		decimal _z;     /* Counts to estimate the future parameter (used by EM). */
		number _cfreq;  /* Number of data samples classified on this component.  */
	}gauss;

	typedef struct{
		gauss *mix;       /* Vector of gaussian distributions (above struct).    */
		number dimension; /* Number of dimensions of the gaussian mixture model. */
		number num;       /* Number of components of the gaussian mixture model. */
		decimal *mcov;    /* Minimum allowed covariances to avoid singularities. */
		decimal llh;      /* LogLikelihood after the training with EM algorithm. */
	}gmm;

	/* Public function prototypes to work with Gaussian Mixture Models. */
	gmm *gmm_load(char*);
	gmm *gmm_initialize(data*,number);
	void gmm_save(char*,gmm*);
	void gmm_save_log(char*,gmm*);
	void gmm_delete(gmm*);
	decimal gmm_classify(char*,data*,gmm*,gmm*,workers*);
	decimal gmm_simple_classify(data*,gmm*,gmm*,workers*);
	decimal gmm_EMtrain(data*,gmm*,workers*);
	gmm *gmm_merge(gmm*,data*,decimal,workers*);
	/* Private function prototypes to work with Gaussian Mixture Models. */
	gmm *gmm_create(number,number);
	void gmm_init_classifier(gmm*);

#endif
