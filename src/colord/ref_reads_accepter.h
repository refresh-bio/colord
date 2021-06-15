/*******************************************************************************
 
 CoLoRd 
 Copyright (C) 2021, M. Kokot, S. Deorowicz, and A. Gudys
 https://github.com/refresh-bio/CoLoRd

 This program is free software: you can redistribute it and/or modify it under 
 the terms of the GNU General Public License as published by the Free Software 
 Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY 
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
 A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this 
 program. If not, see https://www.gnu.org/licenses/.

******************************************************************************/
#pragma once
#include "defs.h"
#include <random>

class CRefReadsAccepter
{
	const uint32_t range;
	const double exponent;

	std::mt19937 mt;
	std::uniform_real_distribution<double> dist;

	uint32_t n_pseudoreads_form_reference;
public:
	explicit CRefReadsAccepter(uint32_t range, double exponent, uint32_t n_pseudoreads_form_reference) :
		range(range),
		exponent(exponent),
		dist(0.0, 1.0),
		n_pseudoreads_form_reference(n_pseudoreads_form_reference)
	{

	}

	uint32_t GetNAccepted(uint32_t n_reads) const
	{
		CRefReadsAccepter cp(range, exponent, n_pseudoreads_form_reference);
		uint32_t res{};
		for (uint32_t i = 0; i < n_reads + n_pseudoreads_form_reference; ++i)
			res += cp.ShouldAddToReference(i);
		return res;
	}

	bool ShouldAddToReference(uint32_t idx)
	{
		if (idx < n_pseudoreads_form_reference)
			return true;
		uint32_t range_no = (idx - n_pseudoreads_form_reference) / range;
		double accept_prob = pow(1.0 / (range_no + 1ul), exponent);
		return dist(mt) <= accept_prob;
	}
};