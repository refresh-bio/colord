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

#include <string>
#include <vector>
#include <cctype>
#include <algorithm>
#include <iostream>
#include <cassert>
#include "edlib.h"


//#define VALIDATE_REFACTOR_EDIT_SCRIPT

//debug helper funciton
//#if 0
inline void printWithAlignment(std::string_view seq, std::string_view editScript)
{
	std::string dispSeq;
	std::string connections;
	std::string result;
	uint32_t input_pos = 0;
	for (size_t i = 0; i < editScript.size(); ++i)
	{
		if (editScript[i] == 'M') //match
		{
			connections.push_back('|');
			dispSeq.push_back(seq[input_pos]);
			result.push_back(seq[input_pos++]);
		}
		else if (editScript[i] == 'D')		//deletion
		{
			dispSeq.push_back(seq[input_pos]);
			++input_pos;
			connections.push_back(' ');
			result.push_back('-');
		}

		//insertion
		else if (editScript[i] == 'A')
			result.push_back('A'), dispSeq.push_back('-'), connections.push_back(' ');
		else if (editScript[i] == 'C')
			result.push_back('C'), dispSeq.push_back('-'), connections.push_back(' ');
		else if (editScript[i] == 'G')
			result.push_back('G'), dispSeq.push_back('-'), connections.push_back(' ');
		else if (editScript[i] == 'T')
			result.push_back('T'), dispSeq.push_back('-'), connections.push_back(' ');

		//substitution
		else
		{
			dispSeq.push_back(seq[input_pos]);
			connections.push_back(' ');

			char symb = CMissmatchCoder::decode_mismatch_symb(seq[input_pos], editScript[i]);
			++input_pos;

			result.push_back(symb);
		}
	}
	std::cerr << dispSeq << "\n";
	std::cerr << connections << "\n";
	std::cerr << result << "\n";
}


inline std::string decode_edit_script(read_view input, std::string_view edit_script, uint32_t& input_pos)
{
	std::string result;
	input_pos = 0;
	for (size_t i = 0; i < edit_script.size(); ++i)
	{
		if (edit_script[i] == 'M') //match
			result.push_back("ACGT"[input[input_pos++]]);
		else if (edit_script[i] == 'D')		//deletion
			++input_pos;
		
		//insertion
		else if(
			edit_script[i] == 'A' || 
			edit_script[i] == 'C' || 
			edit_script[i] == 'G' || 
			edit_script[i] == 'T'
			)
			result.push_back(edit_script[i]);
		
		//substitution
		else
		{
			char symb = CMissmatchCoder::decode_mismatch_symb("ACGT"[input[input_pos]], edit_script[i]);
			++input_pos;			
			result.push_back(symb);
		}
	}
	return result;
}

inline std::string decode_edit_script(read_view input, std::string_view edit_script)
{
	uint32_t input_pos;
	return decode_edit_script(input, edit_script, input_pos);
}

//assert purposes. Check if changes contained in the edit script are in fact of given edit distance
inline bool check_edit_script(std::string_view edit_script, uint32_t edit_distance)
{
	uint32_t ed{};
	for (size_t i = 0; i < edit_script.size(); ++i)
		if (edit_script[i] != 'M')
			++ed;
	return ed == edit_distance;
}
//#endif


class EditDistanceMemoryManager
{
	std::vector<uint32_t> memory;
	uint32_t n, m;
public:
	void SetDimmension(uint32_t _n, uint32_t _m)
	{
		n = _n;
		m = _m;
		if (memory.size() < static_cast<uint64_t>(n) * m)
			memory.resize(static_cast<uint64_t>(n) * m);
	}
	uint32_t* operator[](size_t index)
	{
		return memory.data() + m * index;
	}
};


//returns edit script and edit distance
//currently called only if s1 < 2 || s2 < 2, because it seems that edlib does not hanlde such trivial cases correctly
//This function shuld be rewritten to return 0 if one of input is of length 0, appropriate edit script should be returnet (only insertions and deletions)
//for lenght = 1 find if it occur in the second one and return appropriate result
inline std::pair<std::string, uint32_t> find_edit_dist(read_view s1, read_view s2)
{	
	EditDistanceMemoryManager c;
	c.SetDimmension(static_cast<uint32_t>(s1.length() + 1), static_cast<uint32_t>(s2.length() + 1));

	const uint32_t ci = 1; //insertion and delection cost	
	for (uint32_t i = 0; i < s1.length() + 1; ++i)
		c[i][0] = i * ci;
	for (uint32_t j = 0; j < s2.length() + 1; ++j)
		c[0][j] = j * ci;

	for (uint32_t i = 1; i < s1.length() + 1; ++i)
		for (uint32_t j = 1; j < s2.length() + 1; ++j)
		{
			auto _a = c[i][j - 1] + ci;
			auto _b = c[i - 1][j] + ci;
			auto _c = c[i - 1][j - 1] + (s1[i - 1] != s2[j - 1]);
			c[i][j] = min3(_a, _b, _c);
		}
		
	//backtracking
	uint64_t index_1 = s1.length();
	uint64_t index_2 = s2.length();
	//D - deletion
	//a|c|g|t - insertion of symbol
	//A|C|G|T - substitute to symbol
	//M - match

	std::string conversion_required;

	while (index_1 > 0 && index_2 > 0)
	{
		uint32_t from_up = c[index_1 - 1][index_2] + ci;
		uint32_t from_left = c[index_1][index_2 - 1] + ci;
		uint32_t from_diag = c[index_1 - 1][index_2 - 1] + (s1[static_cast<int32_t>(index_1 - 1)] != s2[static_cast<int32_t>(index_2 - 1)]);
		if (from_up == c[index_1][index_2])
		{
		conversion_required.push_back('D');
		--index_1;
		}
		else if (from_left == c[index_1][index_2])
		{
			conversion_required.push_back("ACGT"[s2[static_cast<int32_t>(index_2 - 1)]]);
			--index_2;
		}
		else if (from_diag == c[index_1][index_2])
		{
			if (s1[static_cast<int32_t>(index_1 - 1)] == s2[static_cast<int32_t>(index_2 - 1)])
			{
				conversion_required.push_back('M');
			}
			else
			{
				conversion_required.push_back(CMissmatchCoder::encode_missmatch_symb("ACGT"[s1[static_cast<int32_t>(index_1 - 1)]], "ACGT"[s2[static_cast<int32_t>(index_2 - 1)]]));
			}
			--index_1;
			--index_2;
		}
		else
		{
			std::cerr << "Error: cirtical, contact authors, file: " << __FILE__ << ", line: " << __LINE__ << "\n";
			exit(1);
		}
	}
	while (index_1 > 0)
	{
		conversion_required.push_back('D');
		--index_1;
	}
	while (index_2 > 0)
	{	
		conversion_required.push_back("ACGT"[s2[static_cast<int32_t>(index_2 - 1)]]);
		--index_2;
	}

	std::reverse(std::begin(conversion_required), std::end(conversion_required));

	assert(check_edit_script(conversion_required, c[s1.size()][s2.size()]));

	//auto dec = decode_edit_script(s1, conversion_required);
	//assert(s2 == dec);

	return { conversion_required, c[s1.size()][s2.size()] };
}



struct EditDistRes
{
	std::string editScript;
	uint32_t editDist;
};


inline EditDistRes get_edit_dist_on_seq_empty(read_view refPart, read_view encPart)
{
	EditDistRes ed;
	if (refPart.size() == 0)
	{
		ed.editDist = static_cast<uint32_t>(encPart.size());
		ed.editScript.reserve(ed.editDist);
		for (auto c : encPart)
			ed.editScript.push_back("ACGT"[c]);
		return ed;
	}
	else //encPart.size() == 0
	{
		ed.editDist = static_cast<uint32_t>(refPart.size());
		ed.editScript.assign(ed.editDist, 'D');
	}
	return ed;
}

/*
* sequences may not be empty!
*/
inline EditDistRes find_edit_dist_with_edlib_ex(read_view ref, read_view enc, EdlibAlignMode mode = EDLIB_MODE_NW)
{
	EditDistRes res;
	if (ref.length() < 2 || enc.length() < 2 || (ref.length() < 15 && enc.length() < 15 && mode == EDLIB_MODE_NW))
	{
		auto tmp = find_edit_dist(ref, enc); 
		res.editDist = tmp.second;
		res.editScript = tmp.first;
		return res;
	}

	EdlibAlignResult result = edlibAlign((const char*)ref.data(), static_cast<int>(ref.length()), (const char*)enc.data(), static_cast<int>(enc.length()),
		edlibNewAlignConfig(-1, mode, EDLIB_TASK_PATH, NULL, 0)
	);
	std::string conversion_required;
	if (result.status == EDLIB_STATUS_OK)
	{
		//0 - match
		//1 - insertion to target, enc
		//2 - insertion to query, ref
		//3 - missmatch		
		uint32_t pos_ref{}, pos_enc{};
		conversion_required.reserve(result.alignmentLength);
		for (int i = 0; i < result.alignmentLength; ++i)
		{
			char c = result.alignment[i];
			if (c == 0)
			{
				++pos_ref, ++pos_enc;
				conversion_required.push_back('M');
			}
			else if (c == 1)
			{
				conversion_required.push_back('D');
				++pos_ref;
			}
			else if (c == 2)				
				conversion_required.push_back("ACGT"[enc[pos_enc++]]);
			else if (c == 3)
			{
				conversion_required.push_back(CMissmatchCoder::encode_missmatch_symb("ACGT"[ref[pos_ref]], "ACGT"[enc[pos_enc]]));
				++pos_ref, ++pos_enc;
			}
		}
		edlibFreeAlignResult(result);

		res.editScript = move(conversion_required);
		res.editDist = result.editDistance;
		return res;

	}
	else
	{
		std::cerr << "edlib error\n";
		exit(1);
	}
}

/*
* sequences may not be empty!
*/
inline EditDistRes find_edit_dist_with_edlib_ex_odwr(read_view ref, read_view enc, uint32_t& ref_end, EdlibAlignMode mode = EDLIB_MODE_NW)
{
	EditDistRes res;
	if (ref.length() < 2 || enc.length() < 2)
	{
		auto tmp = find_edit_dist(ref, enc);
		res.editDist = tmp.second;
		res.editScript = tmp.first;
		ref_end = static_cast<uint32_t>(ref.size() - 1);
		return res;
	}

	EdlibAlignResult result = edlibAlign((const char*)enc.data(), static_cast<int>(enc.length()), (const char*)ref.data(), static_cast<int>(ref.length()),
		edlibNewAlignConfig(-1, mode, EDLIB_TASK_PATH, NULL, 0)
	);
	std::string conversion_required;
	if (result.status == EDLIB_STATUS_OK)
	{
		assert(result.numLocations > 0);
		ref_end = result.endLocations[0];
		//0 - match
		//1 - insertion to target, ref  (deletion in query, enc)
		//2 - insertion to query, enc (deletion in target, ref)
		//3 - missmatch		
		uint32_t pos_ref{}, pos_enc{};
		conversion_required.reserve(result.alignmentLength);
		for (int i = 0; i < result.alignmentLength; ++i)
		{
			char c = result.alignment[i];
			if (c == 0)
			{
				++pos_ref, ++pos_enc;
				conversion_required.push_back('M');
			}
			else if (c == 1)
			{
				conversion_required.push_back("ACGT"[enc[pos_enc++]]);
			}
			else if (c == 2)
			{
				conversion_required.push_back('D');
				++pos_ref;
			}
			else if (c == 3)
			{
				conversion_required.push_back(CMissmatchCoder::encode_missmatch_symb("ACGT"[ref[pos_ref]], "ACGT"[enc[pos_enc]]));
				++pos_ref, ++pos_enc;
			}

		}
		edlibFreeAlignResult(result);

		res.editScript = move(conversion_required);
		res.editDist = result.editDistance;
		return res;

	}
	else
	{
		std::cerr << "edlib error\n";
		exit(1);
	}
}

/*
* sequences may not be empty!
*/
inline EditDistRes find_edit_dist_with_edlib_ex_odwr_reverse(read_view ref, read_view enc, uint32_t max_symb_in_ref, uint32_t& ref_offset, EdlibAlignMode mode = EDLIB_MODE_NW)
{
	read_t _ref = ref.get_reversed();	
	read_t _enc = enc.get_reversed();
	
	uint32_t ref_end;
	auto res = find_edit_dist_with_edlib_ex_odwr(read_view(_ref).substr(0, max_symb_in_ref), _enc, ref_end, mode);
	std::reverse(res.editScript.begin(), res.editScript.end());

	auto max_end = static_cast<uint32_t>(ref.size() - 1);
	ref_offset = max_end - ref_end;

	return res;
}


inline bool TheSameSymbolsInRange(read_view input, uint32_t start, uint32_t end)
{
	for(uint32_t i = start + 1 ; i < end ; ++i)
		if(input[start] != input[i])
			return false;
	return true;
}

inline bool NoMissmatchInRange(const std::string& input, uint32_t start, uint32_t end)
{
	for(uint32_t i = start ; i < end ; ++i)
		if(CMissmatchCoder::is_missmatch_es_symb(input[i]))
			return false;
	return true;
}

inline void FixInRange(std::string& edit_script, uint64_t start, uint64_t end)
{
	if(end < start + 2)
		return;
	--end;

	while(true)
	{
		while(start < end && edit_script[start] == 'M') ++start;
		while(start < end && edit_script[end] != 'M') --end;
		if(start == end)
			break;
		std::swap(edit_script[start], edit_script[end]);
	}
}


#ifdef VALIDATE_REFACTOR_EDIT_SCRIPT

/*
	edit script is refactored OK if deletions are ok and insertions are ok
	example of not ok deletions:
	reference:   TTT
	edit script: DMM
	example of ok deletions
	reference: TTT
	edit script: MMD
	in summary for:
	1. a range of identical symbols in reference and 
	2. corresponding range in edit script containing only deletions and matches 
	deletions should be at the end of edit scirpt range

	for insertion by analogy
*/
bool is_refactored_OK(read_view input, std::string_view edit_script)
{
	uint32_t input_pos = 0;
	for (size_t i = 0; i < edit_script.size()-i; ++i)
	{
		//check if delections refactored ok
		if (edit_script[i] == 'D' && edit_script[i + 1] == 'M' && input[input_pos] == input[input_pos + 1])
			return false;

		//check if insertions refactored ok
		if (edit_script[i + 1] == 'M')
		{
			char s = edit_script[i];
			if (s == 'A' || s == 'C' || s == 'G' || s == 'T')
			{
				if (s == "ACGT"[input[input_pos]])
					return false;
			}
		}

		bool insertion = edit_script[i] == 'A' || edit_script[i] == 'C' || edit_script[i] == 'G' || edit_script[i] == 'T';
		if (!insertion)
			++input_pos;		
	}
	return true;
}


void validate_refactored_edit_scirpt(const std::string& org_es, const std::string& refactored_es, read_view ref_part)
{
	if (org_es == refactored_es)	
		return;
	
	auto dec = decode_edit_script(ref_part, refactored_es);
	
	if (dec != decode_edit_script(ref_part, org_es))
	{
		std::cerr << "Error: refactored edit script fatal error\n";
		exit(1);
	}

	
	if (!is_refactored_OK(ref_part, refactored_es))		
		std::cerr << "Warning: refactor edit script does not work correctly!\n";			
}



bool is_refactored_dels(read_view input, std::string_view edit_script)
{
	uint32_t input_pos = 0;
	for (size_t i = 0; i < edit_script.size() - i; ++i)
	{
		if (edit_script[i] == 'D' && edit_script[i + 1] == 'M' && input[input_pos] == input[input_pos + 1])
		{
			std::cerr << "del error es pos: " << i << "\n";
			std::cerr << "del error ref pos: " << input_pos << "\n";
			return false;
		}

		if (edit_script[i] == 'M') //match
			input_pos++;
		else if (edit_script[i] == 'D')		//deletion
			++input_pos;

		//insertion
		else if (
			edit_script[i] == 'A' ||
			edit_script[i] == 'C' ||
			edit_script[i] == 'G' ||
			edit_script[i] == 'T'
			)
			;

		//substitution
		else
		{
			char symb = CMissmatchCoder::decode_mismatch_symb(input[input_pos], edit_script[i]);
			++input_pos;
		}
	}
	return true;
}


bool is_refactored_ins(read_view input, std::string_view edit_script)
{
	uint32_t input_pos = 0;
	for (size_t i = 0; i < edit_script.size() - i; ++i)
	{
		if (edit_script[i + 1] == 'M')
		{
			char s = edit_script[i];
			if (s == 'A' || s == 'C' || s == 'G' || s == 'T')
			{
				if (s == "ACGT"[input[input_pos]])
					return false;
			}
		}
		if (edit_script[i] == 'M') //match
			input_pos++;
		else if (edit_script[i] == 'D')		//deletion
			++input_pos;

		//insertion
		else if (
			edit_script[i] == 'A' ||
			edit_script[i] == 'C' ||
			edit_script[i] == 'G' ||
			edit_script[i] == 'T'
			)
			;

		//substitution
		else
		{
			char symb = CMissmatchCoder::decode_mismatch_symb(input[input_pos], edit_script[i]);
			++input_pos;
		}
	}
	return true;
}

#endif

inline void refactor_edit_script_dels(read_view ref, read_view enc, std::string& edit_script)
{
	uint32_t ref_start = 0;
	uint32_t ref_pos = 0;
	uint32_t es_start = 0;
	for (uint32_t es_pos = 0; es_pos < edit_script.length(); ++es_pos)
	{
		auto es_s = edit_script[es_pos];

//		bool match = es_s == 'M';
//		bool del = es_s == 'D';		
		bool missmatch = CMissmatchCoder::is_missmatch_es_symb(es_s);
		bool ins = es_s == 'A' || es_s == 'C' || es_s == 'G' || es_s == 'T';

//		assert(match || del || missmatch || ins);

		if (ins || missmatch || ref[ref_start] != ref[ref_pos])
		{
			FixInRange(edit_script, es_start, es_pos);

			es_start = es_pos;
			if (ins || missmatch)
				++es_start;
			ref_start = ref_pos;
		}

		if (!ins)
			++ref_pos;
	}

	FixInRange(edit_script, es_start, edit_script.length());
}

inline void refactor_edit_script_ins(read_view ref, read_view enc, std::string& edit_script)
{
//	auto es_copy = edit_script;
	uint32_t enc_start = 0;
	uint32_t enc_pos = 0;
	uint32_t es_start = 0;
	for (uint32_t es_pos = 0; es_pos < edit_script.length(); ++es_pos)
	{
		auto es_s = edit_script[es_pos];

//		bool match = es_s == 'M';
		bool del = es_s == 'D';		
		bool missmatch = CMissmatchCoder::is_missmatch_es_symb(es_s);
//		bool ins = es_s == 'A' || es_s == 'C' || es_s == 'G' || es_s == 'T';

//		assert(match || del || missmatch || ins);

		if (del || missmatch || enc[enc_start] != enc[enc_pos])
		{
			FixInRange(edit_script, es_start, es_pos);

			es_start = es_pos;
			if (del || missmatch)
				++es_start;
			enc_start = enc_pos;
		}

		if (!del)
			++enc_pos;
	}

	FixInRange(edit_script, es_start, edit_script.length());
}




inline void refactor_edit_script(read_view ref, read_view enc, std::string& edit_script)
{
#ifdef VALIDATE_REFACTOR_EDIT_SCRIPT
	std::string es_copy = edit_script;
#endif 
	refactor_edit_script_dels(ref, enc, edit_script);
	refactor_edit_script_ins(ref, enc, edit_script);
#ifdef VALIDATE_REFACTOR_EDIT_SCRIPT
	validate_refactored_edit_scirpt(es_copy, edit_script, ref);
#endif
}
