/*先頭で日本語を打ち込んでおけばソースツリーで表示したときに文字化けしないらしいので*/

#ifndef __INFO__
#define __INFO__

#include "common.h"
#include "nariinfocontroller.h"
#include "narifile.h"
#include "naricommon.h"
#include <string>

struct info
{
	std::string dir_Ref;
	std::string dir_mvd;
	std::string dir_Ans;
	std::string dir_out;
	std::string dir_list;
	std::string case_flist;
	std::string case_rlist;
	int tmp;
	int rangex;
	int rangey;
	int rangez;

	inline void input(const std::string &path)
	{
		nari::infocontroller info;
		info.load(path);
		dir_Ref = nari::file::add_delim(info.get_as_str("dir_Ref"));
		dir_mvd = nari::file::add_delim(info.get_as_str("dir_mvd"));
		dir_Ans = nari::file::add_delim(info.get_as_str("dir_Ans"));
		dir_out = nari::file::add_delim(info.get_as_str("dir_out"));
		dir_list = nari::file::add_delim(info.get_as_str("dir_txt"));
		case_flist = info.get_as_str("case_f");
		case_rlist = info.get_as_str("case_r");
		tmp = info.get_as_int("tmp_size");
		rangex = info.get_as_int("range_x");
		rangey = info.get_as_int("range_y");
		rangez = info.get_as_int("range_z");

	}
};


#endif