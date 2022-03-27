	for (k = 1; k < elem_d_ps_nsoil[tid]; k++)  // 不确定的循环次数
	{
		hcpct = elem_d_ws_sh2o[k][tid] * CH2O + (1.0 - elem_d_soil_smcmax[tid]) * csoil_loc +
			(elem_d_soil_smcmax[tid] - elem_d_ws_smc[k][tid]) * CP +
			(elem_d_ws_smc[k][tid] - elem_d_ws_sh2o[k][tid]) * CPICE;

		/* This section for Layer 2 or greater, but not last layer.
		* Calculate thermal diffusivity for this layer. */
		if (k != elem_d_ps_nsoil[tid] - 1)
		{
			/* Calc the vertical soil temp gradient thru this layer */
			df1n = TDfCnd(elem_d_ws_smc[k][tid], elem_d_soil_quartz[tid], elem_d_soil_smcmax[tid], 
				elem_d_soil_smcmin[tid],
				elem_d_ws_sh2o[k][tid]);

			/* Urban */
			df1n = (elem_d_lc_isurban[tid]) ? 3.24 : df1n;

			denom = 0.5 * (elem_d_ps_zsoil[k - 1][tid] - elem_d_ps_zsoil[k + 1][tid]);

			/* Calc the matrix coef, ci, after calc'ng its partial product */
			dtsdz2 = (elem_d_es_stc[k][tid] - elem_d_es_stc[k + 1][tid]) / denom;
			ddz2 = 2.0 / (elem_d_ps_zsoil[k - 1][tid] - elem_d_ps_zsoil[k + 1][tid]);

			/* If temperature averaging invoked (itavg = true; else skip)
			* Calculate temp at bottom of layer. */
			ci[k] = -df1n * ddz2 / ((elem_d_ps_zsoil[k - 1][tid] - elem_d_ps_zsoil[k][tid]) * hcpct);
			if (itavg)
			{
				tbk1 = TBnd(elem_d_es_stc[k][tid], elem_d_es_stc[k + 1][tid], 
					elem_d_ps_zsoil[tid], elem_d_ps_zbot[tid], k,
					elem_d_ps_nsoil[tid]);

			}
		}  // if
		else
		{
			/* Special case of bottom soil layer
			* Calculate thermal diffusivity for bottom layer. */
			df1n = TDfCnd(elem_d_ws_smc[k][tid], elem_d_soil_quartz[tid], elem_d_soil_smcmax[tid], 
				elem_d_soil_smcmin[tid],
				elem_d_ws_sh2o[k][tid]);

			/* Urban */
			df1n = (elem_d_lc_isurban[tid]) ? 3.24 : df1n;

			/* Calc the vertical soil temp gradient thru bottom layer. */
			denom = 0.5 * (elem_d_ps_zsoil[k - 1][tid] + elem_d_ps_zsoil[k][tid]) - elem_d_ps_zbot[tid];
			dtsdz2 = (elem_d_es_stc[k][tid] - elem_d_ps_tbot[tid]) / denom;

			/* Set matrix coef, ci to zero if bottom layer. */
			ci[k] = 0.0;

			/* If temperature averaging invoked (itavg = true; else skip)
			* Calculate temp at bottom of last layer. */
			if (itavg)
			{
				tbk1 = TBnd(elem_d_es_stc[k][tid], elem_d_ps_tbot[tid], 
					elem_d_ps_zsoil[tid], elem_d_ps_zbot[tid], k,
					elem_d_ps_nsoil[tid]);

			}    /* This ends special loop for bottom layer. */
		}  // if( k != elem_d_ps_nsoil[tid] - 1)

		/* Calculate rhsts for this layer after calc'ng a partial product. */
		denom = (elem_d_ps_zsoil[k][tid] - elem_d_ps_zsoil[k - 1][tid]) * hcpct;
		rhsts[k] = (df1n * dtsdz2 - df1k * dtsdz) / denom;
		qtot = -1.0 * denom * rhsts[k];

		sice = elem_d_ws_smc[k][tid] - elem_d_ws_sh2o[k][tid];

		if (itavg)
		{
		//	tavg = TmpAvg(tbk, elem_d_es_stc[k][tid], tbk1, elem_d_ps_zsoil[tid], k);
			tavg = TmpAvg(tbk, elem_d_es_stc[k][tid], tbk1, k,
			    elem_d_ps_zsoil0[tid], elem_d_ps_zsoil1[tid], elem_d_ps_zsoil2[tid],
				elem_d_ps_zsoil3[tid], elem_d_ps_zsoil4[tid], elem_d_ps_zsoil5[tid],
				elem_d_ps_zsoil6[tid], elem_d_ps_zsoil7[tid], elem_d_ps_zsoil8[tid],
				elem_d_ps_zsoil9[tid], elem_d_ps_zsoil10[tid]);

			if ((sice > 0.0) || (elem_d_es_stc[k][tid] < TFREEZ) || (tbk < TFREEZ) ||(tbk1 < TFREEZ))	
			{
				SnkSrc(tid, &tsnsr, tavg, elem_d_ws_smc[k][tid], &elem_d_ws_sh2o[k][tid],
					elem_d_ps_zsoil[tid],
					dt, k, qtot,
					elem_d_soil_alpha,
					elem_d_soil_beta,
					elem_d_soil_smcmin,
					elem_d_soil_smcmax
					);
				rhsts[k] = rhsts[k] - tsnsr / denom;
			}
		}
		else
		{
			if ((sice > 0.0) || (elem_d_es_stc[k][tid] < TFREEZ))
			{
				SnkSrc(tid, &tsnsr, elem_d_es_stc[k][tid], elem_d_ws_smc[k][tid],
					&elem_d_ws_sh2o[k][tid],
					elem_d_ps_zsoil[tid], dt, k, qtot,
					elem_d_soil_alpha,
					elem_d_soil_beta,
					elem_d_soil_smcmin,
					elem_d_soil_smcmax
					);
				rhsts[k] = rhsts[k] - tsnsr / denom;
			}
		}

		/* Calc matrix coefs, ai, and bi for this layer. */
		ai[k] = -df1k * ddz / ((elem_d_ps_zsoil[k - 1][tid] - elem_d_ps_zsoil[k][tid]) * hcpct);
		bi[k] = -(ai[k] + ci[k]);

		/* Reset values of df1, dtsdz, ddz, and tbk for loop to next soil layer.
		*/
		tbk = tbk1;
		df1k = df1n;
		dtsdz = dtsdz2;
		ddz = ddz2;
	}  // for (k=1; )