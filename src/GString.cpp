
//
//  GString.cpp
//  GFC
//
//  Created by lizhen on 15/9/21.
//  Copyright © 2015年 lizhen. All rights reserved.
//

//============================================================================
//
//  This file is part of GFC, the GNSS FOUNDATION CLASS.
//
//  The GFC is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 3.0 of the License, or
//  any later version.
//
//  The GFC is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with GFC; if not, write to the Free Software Foundation,
//  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110, USA
//
//  Copyright 2015, lizhen
//
//============================================================================


#include "GString.h"

namespace gfc
{
    
    void GString::hexDumpData(std::ostream& s,
                              const GString& data,
                              unsigned indent,
                              HexDumpDataConfig cfg)
    {
        GString instr(indent, ' ');
        hexDumpData(s, data, instr, cfg);
    }
    
    void GString::hexDumpData(std::ostream& s, const GString& data,
                              const GString& tag, HexDumpDataConfig cfg)
    {
        GString ascii="";
        
        uint_t indent =  static_cast<uint_t>(tag.length()); //(uint_t)tag.length();
        int_t col = 0;
        
        int_t datasize = static_cast<uint_t>(data.size());
        
        GString groupws(cfg.groupWS, ' ');
        GString  group2ws(cfg.group2WS, ' ');
        GString  indexws(cfg.indexWS, ' ');
        GString  textws(cfg.textWS, ' ');
        uint_t linesize;
        
        if (cfg.groupBy && ((cfg.bytesPerLine % cfg.groupBy) != 0))
        {
            s << "hexDumpData: cfg.bytesPerLine % cfg.groupBy != 0"
            << std::endl;
            return;
        }
        if (cfg.group2By && ((cfg.bytesPerLine % cfg.group2By) != 0))
        {
            s << "hexDumpData: cfg.bytesPerLine % cfg.group2By != 0"
            << std::endl;
            return;
        }
        if (cfg.groupBy && ((cfg.group2By % cfg.groupBy) != 0))
        {
            s << "hexDumpData: cfg.group2By % cfg.groupBy != 0"
            << std::endl;
            return;
        }
        
        // line format:
        // <tag><index>:<indexws><group1byte1>...<group1byte[groupBy]><groupws>...<group[group2By]byte1>...<group[group2By]byte[groupBy]><group2ws>....<byte[bytesPerLine]><textws><separator><text><separator>\n
        linesize = indent;
        if (cfg.showIndex)
            linesize += cfg.idxDigits + 1 + cfg.indexWS;
        linesize += cfg.bytesPerLine * 2;
        unsigned w2 = 0;
        unsigned w1 = 0;
        if (cfg.group2By)
            w2 = (cfg.bytesPerLine / cfg.group2By) - 1;
        if (cfg.groupBy)
            w1 = (cfg.bytesPerLine / cfg.groupBy) - w2 - 1;
        if (cfg.groupBy > 0)
            linesize += cfg.groupWS * w1;
        if (cfg.group2By > 0)
            linesize += cfg.group2WS * w2;
        /*
         linesize doesn't include text stuff
         if (cfg.showText)
         linesize += cfg.textWS + cfg.bytesPerLine;
         if (cfg.separator)
         linesize += 2;
         */
        
        for (int_t i=0; i<datasize; i++)
        {
            if (i%cfg.bytesPerLine==0)
            {
                s << tag;
                col = indent;
                if (cfg.showIndex)
                {
                    if (cfg.hexIndex)
                    {
                        s << std::hex;
                        if (cfg.upperHex)
                            s << std::uppercase;
                        else
                            s << std::nouppercase;
                    }
                    else
                        s << std::dec;
                    s << std::setfill('0');
                    s << std::setw(cfg.idxDigits) << i << ":" << indexws;
                    s << std::dec << std::nouppercase;
                }
                col += cfg.idxDigits + 1 + cfg.indexWS;
            }
            
            unsigned char c = data[i];
            
            if (isprint(c))
                ascii += c;
            else
                ascii += '.';
            if (cfg.upperHex)
                s << std::uppercase;
            else
                s << std::nouppercase;
            s << std::hex << std::setw(2) << (int_t)c << std::dec
            << std::nouppercase;
            col += 2;
            if (((i % cfg.bytesPerLine) == (cfg.bytesPerLine-1)) ||
                (i == (datasize-1)))
            {
                if (cfg.showText)
                {
                    int_t extra = linesize-col;
                    GString space(extra, ' ');
                    s << space << textws;
                    if (cfg.separator)
                        s << cfg.separator;
                    s << ascii;
                    if (cfg.separator)
                        s << cfg.separator;
                    s << std::endl;
                }
                // this *should* be updated at the beginning of the loop
                //col=indent+6;
                ascii.erase();
            }
            else if (cfg.group2By && ((i % cfg.group2By) == (cfg.group2By-1)))
            {
                s << group2ws;
                col += cfg.group2WS;
            }
            else if (cfg.groupBy && ((i % cfg.groupBy) == (cfg.groupBy-1)))
            {
                s << groupws;
                col += cfg.groupWS;
            }
        }
    }
    
    
    // Keep searching for aString at the start of s
    // until num == 0 or aString is not found at the start of s
    ////从头开始去掉最开始的字符串s
    //num 用于处理重复的次数，默认是去掉所有重复的aString
    GString GString::stripLeading(const GString& aString, GString::size_type num)
    throw(GStringException)
    {
        try
        {
            GString s = *this;
            if (aString == "")
                return s;
            while( (num > 0)
                  &&(  s.find(aString,0) == 0)
                  &&(  s.length() > 0))
            {
                s.erase(0,aString.length());
                --num;
            }
            return s;
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    //返回值为空，表示对字符串本身进行处理，处理结果将当前字符串的结果覆盖掉
    GString GString::stripLeading_v(const GString& aString, GString::size_type num)
    throw(GStringException)
    {
        try
        {
            if (aString == "")
                return *this;
            while( (num > 0)
                  &&(  this->find(aString,0) == 0)
                  &&(  this->length() > 0))
            {
                this->erase(0,aString.length());
                --num;
            }
            return *this;
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    
    GString GString::stripTrailing( const GString& aString, GString::size_type num )
    throw(GStringException)
    {
        try
        {
            GString s = *this;
            GString::size_type pos = this->length() - aString.length();
            
            // empty string, etc.
            if ((pos > s.length()) || (aString == ""))
                return s;
            
            while((num > 0) &&
                  (s.rfind(aString,pos) == pos) &&
                  (s.length() > 0))
            {
                s.erase(pos, GString::npos);
                --num;
                pos = s.length() - aString.length();
            }
            return s;
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + std::string(e.what()));
            GFC_THROW(strexc);
        }
        
    }
    
    GString GString::stripTrailing_v( const GString& aString, GString::size_type num )
    throw(GStringException)
    {
        try
        {
            GString::size_type pos = this->length() - aString.length();
            
            // empty string, etc.
            if ((pos > this->length()) || (aString == ""))
                return *this;
            
            while((num > 0) &&
                  (this->rfind(aString,pos) == pos) &&
                  (this->length() > 0))
            {
                this->erase(pos, GString::npos);
                --num;
                pos = this->length() - aString.length();
            }
            return *this;
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + std::string(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    GString GString::strip(  const GString& aString, GString::size_type num)
    throw(GStringException)
    {
        //不修改本身字符串就得使用中间变量
        GString s1 = stripLeading(aString, num);
        GString s2 = s1.stripTrailing(aString, num);
        return s2;
    }
    
    GString GString::strip_v(  const GString& aString, GString::size_type num)
    throw(GStringException)
    {
        //连续修改本身字符串
        stripLeading_v(aString, num);
        stripTrailing_v(aString, num);
        return *this;
    }
    
    GString GString::translate( const GString& inputChars,
                               const GString& outputChars,
                               const char pad)
    {
        GString rv = *this;
        GString::size_type aspos = 0;
        GString::size_type inpos = GString::npos;
        char toc = pad;
        
        // By starting at the last position, we avoid infinite
        // loops in case someone did something dumb, like, for
        // example, setting inputChars=outputChars.
        while ( ( aspos = rv.find_first_of(inputChars, aspos))
               != GString::npos)
        {
            // figure out which char we found;
            inpos = inputChars.find(rv[aspos]);
            if ( outputChars.length()-1 < inpos )
                toc = pad;
            else
                toc = outputChars[inpos];
            rv[aspos] = toc;
            
            aspos++; // try to guarantee no infinite loops
        }
        
        return rv;
    }
    
    GString GString::Greplace(  const GString& originString,
                              const GString& replaceString,
                              GString::size_type startPos ,
                              unsigned numChanges  )
    {
        unsigned count = 0;
        GString::size_type opos = startPos;
        while ( count < numChanges)
        {
            GString::size_type pos = this->find(originString, opos);
            if ( pos != GString::npos )
            {
                count++;
                this->replace(pos, originString.length(), replaceString);
                opos = pos + replaceString.length();
            }
            else
                break;
        }
        return *this;
    }
    
    
    GString GString::rightJustify(   const GString::size_type length, const char pad)
    throw(GStringException)
    {
        try
        {
            GString s = *this;  //复制一份内存，保证不改变this的值
            if( length < s.length())
            {
                s = (s.substr(s.length()-length, GString::npos));  //需要实现赋值重载函数
            }
            else
            {
                s.insert((GString::size_type)0, length-s.length(), pad);
            }
            return s;
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    GString GString::leftJustify(    const GString::size_type length,const char pad )   throw(GStringException)
    {
        try
        {
            GString s = *this;  //复制一份内存，保证不改变this的值
            if( length < s.length() )
            {
                s = s.substr(0, length);
            }
            else
            {
                s.append(length-s.length(), pad);
            }
            return s;
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + std::string(e.what()));
            GFC_THROW(strexc);
        }
        
    }
    
    /**
     * Change the length of a string by adding to the beginning and end.
     * The string \a s is modified to the specified
     * length.  If the string is shorter than
     * \a length, then the string is truncated with the
     * left-most \a length characters remaining.
     * Otherwise, characters are added to the beginning and end of the
     * string until the string is the specified length, where the
     * number of characters added to the beginning and the end
     * does not differ by more than one so the original string
     * is centered.
     * @param s string to be modified.
     * @param length new desired length of string.
     * @param pad character to pad string with (blank by default).
     * @throws StringException if there's a std::exception thrown.
     * @return a reference to \a s.
     */
    GString GString::center( const std::string::size_type length, const char pad) throw(GStringException)
    {
        try
        {
            GString s = *this;
            if( length < s.length())
            {
                s.leftJustify(length,pad);
            }
            else
            {
                GString::size_type leftOff = s.length() + (length - s.length()) / 2;
                s.leftJustify( leftOff, pad);
                s.rightJustify( length, pad);
            }
            return s;
        }
        catch(GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
        
    }
    
    float GString::asFLOAT() throw(GStringException)
    {
        try
        {
            std::string s(this->data());
            std::istringstream is(s);
            float f;
            is >> f;
            return f;
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    long double GString::asLDOUBLE() throw(GStringException)
    {
        try
        {
            std::string s(this->data());
            std::istringstream is(s);
            long double x;
            is >> x;
            return x;
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    /*整型字符串转换为十六进制数据*/
    GString GString::d2x() throw(GStringException)
    {
        try
        {
            //这里需要进行类型检验，必须是整型才能转换
            GString s = *this;
            // remove the integer from s, including
            // leading spaces and 0's
            long l = s.asINT(); //转换为整数
            s.stripLeading();    //去掉空格
            s.stripLeading("0");  //去掉前面的0
            //s.stripLeading(asString<long>(l));
            
            // put the int in a stringstream to convert it
            std::ostringstream st;
            st << std::hex << l << std::dec;
            
            GString tmpstring(st.str());
            return tmpstring.upperCase();
        }
        catch(GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    
    GString GString::x2d()  throw(GStringException)
    {
        try
        {
            GString s = *this;
            // remove the "0x" part, leading zeros and spaces from the
            // string
            // ex. ' 0x003' -> '3'
            s.stripLeading();
            s.stripLeading( "0x", 1);
            s.stripLeading( "0");
            
            // make the stringstream, get the integer, and
            // remove it from the string
            std::istringstream strstr(s);
            int_t i = 0;
            strstr >> std::hex >> i;
            GString aString = s.asString(s.asINT());
            s = s.stripLeading(aString,1);
            //stripLeading(s, asString<int_t>(asInt(s)), 1);
            GString res(i);
            // append the decimal to the existing string
            //s.insert(0,asString<int_t>(i));
            return res;
        }
        catch(GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + std::string(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    
    
    /*
     * 将字符串变为大写字母
     */
    GString  GString::upperCase()
    {
        GString s = *this;
        for( GString::size_type i = 0; i < s.length(); i++)
        {
            s[i] = toupper(s[i]);
        }
        return s;
    }
    
    /*
     *将字符串变为小写字母
     */
    GString GString::lowerCase()
    {
        GString s = *this;
        for( GString::size_type i = 0; i < s.length(); i++)
        {
            s[i] = tolower(s[i]);
        }
        return s;
    }
    
    GString GString::c2x() throw(GStringException)
    {
        const char hexDigits[] = "0123456789ABCDEF";
        try
        {
            GString s = *this;
            std::string old(s);
            const unsigned char *pSource = (unsigned char *)old.c_str();
            unsigned n = static_cast<uint_t>(old.length());
            
            s.resize(n * 2, 0);
            
            for (int_t i = 0; i < (int_t)n * 2;)
            {
                unsigned char c = *pSource++;
                s[i++] = hexDigits[ c / 16 ];
                s[i++] = hexDigits[ c % 16 ];
            }
            
            return s;
        }
        catch( GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch( std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    /*判断一个字符串是否是数字*/
    bool GString::isDigitString()
    {
        if ( size() == 0)
            return false;
        
        const_iterator i = begin();
        
        if( ( *i == '-') || ( *i == '+') )
            i++;
        for( ; i != end(); i++)
            if ( !isdigit(*i) )
                return false;
        return true;
    }
    
    /*判断一个字符串是否是十进制数据的字符串*/
    bool GString::isDecimalString()
    {
        if ( size() == 0)
            return false;
        const_iterator i = begin();
        bool sawdot = false;
        if(( *i == '-') || ( *i == '+'))
            i++;
        for( ; i != end(); i++)
        {
            if ( *i == '.')
            {
                if (sawdot)
                    return false;
                else sawdot = true;
            }
            else if (!isdigit(*i))
                return false;
        }
        return true;
    }
    
    /*判断一个字符串是否是用科学计数法表示的*/
    bool GString::isScientificString()
    {
        if( size() == 0)
            return false;
        
        GString::size_type pos = find_first_of("EeDd");
        if(pos == GString::npos)
            return isDecimalString();
        
        GString mant=substr(0,pos);
        GString exp = substr(pos+1);
        
        return (mant.isDecimalString() && ( exp.size()==0 || exp.isDigitString()));
    }
    
    //判断一个字符串是否全是字母组成
    bool GString::isAlphaString()
    {
        if ( size() == 0)
            return false;
        
        const_iterator i = begin();
        //std::string::size_type index;
        for( ; i != end(); i++)
            if (!isalpha(*i))
                return false;
        return true;
    }
    
    
    GString GString::matches(   const GString& aPattern,
                             const char zeroOrMore,
                             const char oneOrMore,
                             const char anyChar)
		  throw(GStringException)
    {
        GString s = *this;
        GString thisPattern(aPattern);
        GString thisStr(s);
        
        // check if something other than the regex standard
        // characters (*,+,.) is used for those variables
        if ( zeroOrMore != '*')
        {
            //replaceAll(thisPattern, "*", "\\*");
            //replaceAll(thisPattern, GString(1, zeroOrMore), "*");
            thisPattern.Greplace("*","\\*");
            thisPattern.Greplace(GString(1, zeroOrMore), "*");
        }
        if (oneOrMore != '+')
        {
            //replaceAll(thisPattern, "+", "\\+");
            //replaceAll(thisPattern, GString(1, oneOrMore), "+");
            thisPattern.Greplace("+","\\+");
            thisPattern.Greplace(GString(1, oneOrMore), "+");
        }
        if (anyChar != '.')
        {
            //replaceAll(thisPattern, ".", "\\.");
            //replaceAll(thisPattern, GString(1, anyChar), ".");
            thisPattern.Greplace(".","\\.");
            thisPattern.Greplace(GString(1, anyChar), ".");
        }
        
#if defined(_WIN32) && _MSC_VER >= 1700
        try
        {
            std::regex reg (thisPattern);
            
            std::smatch sm;
            if(std::regex_search(thisStr,sm,reg,
                                 std::regex_constants::match_not_bol|
                                 std::regex_constants::match_not_eol))
            {
                return sm.str();
            }
            else
            {
                return std::string();
            }
        }
        catch(std::regex_error& e)
        {
            Exception E(GString("std::regex_error: ") + e.what() );
            GFC_THROW(E);
        }
        
#else
        
        const GString::size_type regErrorBufSize = 512;
        
        regmatch_t matches;
        regex_t regExp;
        char errorMsg[regErrorBufSize];
        int_t rc = regcomp(&regExp, thisPattern.c_str(), REG_EXTENDED);
        
        if (rc != 0)
        {
            regerror(rc, NULL, errorMsg, regErrorBufSize - 1);
            regfree(&regExp);
            GStringException strerr("Regexp error: " + GString(errorMsg));
            GFC_THROW(strerr);
        }
        rc = regexec(&regExp, thisStr.c_str(), 1, &matches, REG_NOTBOL | REG_NOTEOL);
        if ( (rc != 0) && (rc != REG_NOMATCH) )
        {
            regerror(rc, &regExp, errorMsg, regErrorBufSize - 1);
            regfree(&regExp);
            GStringException strerr("Regexp error: " + GString(errorMsg));
            GFC_THROW(strerr);
        }
        
        regfree(&regExp);
        if (rc == REG_NOMATCH)
            return GString();
        else
            return thisStr.substr(matches.rm_so, matches.rm_eo - matches.rm_so);
#endif
    }
    
    template <class T>
    GString GString::formattedPrint( GString& fmt,  GString& pat,  GString& rep, T to ) throw(GStringException)
    {
        
#if defined(_WIN32) && _MSC_VER >= 1700
        
        GString rv(fmt);
        try
        {
            std::regex reg(pat);
            std::smatch m;
            while ( std::regex_search (rv,m,reg) )
            {
                GString mac = m.str();
                mac = mac.Greplace( rep.substr(0,1), rep.substr(1));
                char buffer[1024] = {0};
                sprintf(buffer, mac.c_str(), to);
                rv.replace(m.position(), m.length(), GString(buffer));
            }
        }
        catch(std::regex_error& e)
        {
            Exception E( GString("std::regex_error:")+e.what());
            GFC_THROW(E);
        }
        
        return rv;
#else
        regex_t re;
        const size_t bufferSize = 513;
        char buffer[bufferSize] = {0};
        int rc = regcomp(&re, pat.c_str(), REG_EXTENDED);
        // if the regex doesnt compile, toast =)
        if ( rc != 0)
        {
            regerror(rc, NULL, buffer, bufferSize - 1);
            regfree(&re);
            GStringException se("Regexp error: " + GString(buffer));
            GFC_THROW(se);
        }
        
        regmatch_t r;
        GString rv = fmt;
        while ( regexec(&re, rv.c_str(), 1, &r, 0) == 0 )
        {
            size_t len = r.rm_eo - r.rm_so;
            GString mac = rv.substr(r.rm_so, len);
            mac = mac.Greplace(rep.substr(0,1), rep.substr(1));
            sprintf(buffer, mac.c_str(), to);
            rv.replace(r.rm_so, len, GString(buffer));
        }
        
        regfree(&re);
        return rv;
#endif
    }
    
    GString GString::firstWord( const char delimiter)
    throw(GStringException)
    {
        try
        {
            // return s if there are no delimiters
            GString::size_type pos = find_first_not_of(delimiter);
            if ( pos == GString::npos)
            {
                return *this;
            }
            // find the end delimiter (if any) and return the string
            GString::size_type endPos = find(delimiter, pos);
            if ( endPos == GString::npos)
            {
                return GString(*this,pos, endPos);
            }
            else
            {
                return GString(*this, pos, endPos - pos);
            }
        }
        catch(GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    
    int GString::wordsNum(const char delimiter) throw(GStringException)
    {
        try
        {
            GString t(*this);
            t.stripTrailing_v(delimiter);
            int words = 0;
            while( t.length() )
            {
                t.stripLeading_v( delimiter);
                t.stripLeading_v( t.firstWord(delimiter));
                words++;
            }
            return words;
        }
        catch(GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    
    /*
     *
     *返回从第firstWord往后numWords的字符串，本质上是当前对象的一个字串
     */
    GString GString::wordsString( const GString::size_type firstWord,const GString::size_type numWords,const char delimiter)
    throw(GStringException)
    {
        try
        {
            if ((firstWord == 0) && (numWords == 1))
                return this->firstWord(delimiter);
            if (numWords == 0)
                return "";
            GString::size_type wordNum = 0;
            GString::size_type pos = 0, startPos = 0;
            GString toReturn;
            
            // get position of word wordNum
            pos = find_first_not_of(delimiter, pos);
            while ((pos != GString::npos) && ( pos <= length()))
            {
                if (wordNum == firstWord)
                    startPos = pos;
                // get first delimter after word wordNum
                pos = find(delimiter, pos);
                if (((int)numWords != -1) && ((int)wordNum == (int)(firstWord + (numWords-1))))
                    break;
                pos = find_first_not_of(delimiter, pos);
                wordNum++;
            }
            
            if (pos == GString::npos)  return  substr(startPos);
            
            return  substr(startPos, pos-startPos);
        }
        catch( GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch( std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    /*
     去掉字符串的第一个word
     返去掉后的子串
     本身不发生变化
     */
    GString GString::stripFirstWord(const char delimiter) throw(GStringException)
    {
        try
        {
            GString s = *this;
            s.stripLeading_v(delimiter);
            GString toReturn = s.firstWord(delimiter);
            s.stripLeading_v(toReturn);
            s.stripLeading_v(delimiter);
            return s;
        }
        catch( GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    /*去掉第一个word；本身也会改变*/
    GString GString::stripFirstWord_v(const char delimiter) throw(GStringException)
    {
        try
        {
            stripLeading_v(delimiter);
            GString toReturn = firstWord(delimiter);
            stripLeading_v(toReturn);
            stripLeading_v(delimiter);
            return toReturn;
        }
        catch( GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    std::vector<GString> GString::split(const char delimiter) throw(GStringException)
    {
        try {
            std::vector<GString> rvec;   // vector to return
            GString tempStr(*this);        // copy the input string
            tempStr.stripLeading_v(delimiter); // remove leading delimiters
            while(tempStr.size() > 0)
            {
                rvec.push_back(tempStr.stripFirstWord_v(delimiter));
            }
            return rvec;
        }
        catch(GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    
    GString GString::removeWords( const GString::size_type first,const GString::size_type wordsToReplace,const char delimiter)
    throw(GStringException)
    {
        try
        {
            GString s(*this);
            GString temp(s);
            GString::size_type thisWord;
            // empty out s.  add the new parts of s as they are parsed
            s.erase(0, GString::npos);
            
            // copy the part of the string through word 'first'
            // by appending any delimiters then appending
            // a word for however many words we're keeping.
            for( thisWord = 0; thisWord < first; thisWord++ )
            {
                s.append(temp.find_first_not_of(delimiter),delimiter);
                temp.stripLeading_v(delimiter); //这里需要改变temp的值
                s.append( temp.firstWord() );
                temp.stripLeading_v(temp.firstWord());  //再次修改temp的值
            }
            
            // skip over the number of words to replace, making
            // sure to stop when there's no more string left
            // to skip
            for(thisWord = 0;
                (thisWord < wordsToReplace) &&
                (temp.length() != 0);
                thisWord++)
            {
                temp.stripLeading_v( delimiter);
                temp.stripLeading_v( temp.firstWord() );
            }
            
            // add on any extra words at the end
            s.append(temp);
            
            //决定是否改变本来的值
            *this = s ;
            
            return s;
        }
        catch( GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch( std::exception &e )
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    
    GString GString::double2Sci(  const GString::size_type length,
                                const GString::size_type precision,
                                const GString::size_type explen,
                                bool showPlus)
    {
        // get final exp length, precision and total length
        GString::size_type elen = (explen > 0 ? (explen < 3 ? explen : 3) : 1);
        GString::size_type prec = (precision > 0 ? precision : 1);
        GString::size_type leng = (length > 0 ? length : 1);
        
        // i will be minimum length required with prec==1: force leng if necessary
        int i = (int(leng) - int(elen) - 4);
        if(showPlus) i--;
        if(i > 0 && leng < i) leng = GString::size_type(i);
        
        // set up the stream for writing
        std::stringstream ss;
        ss << std::scientific << std::setprecision((int)prec);
        if(showPlus) ss << std::showpos;
        
        // write GString to a stringstream with precision, sign and in scientific notation
        ss << *this;
        
        // now read that string
        GString str1,str2;
        ss >> str1;
        GString::size_type pos = str1.find_first_of("EDed");    // find exponent
        str2 = str1.substr(0,pos+2);        // str2 = +123.2345e+
        str1 = str1.substr(pos+2);          // str1 = exponent only
        
        // make the exponent length elen
        GString str3(str1.asINT());
        str2+= str3.rightJustify(elen,'0');
        //str2 += rightJustify(StringUtils::asString(StringUtils::asInt(str1)),elen,'0');
        
        // pad if necessary
        if( str2.length() < leng)
        {
            str2 = str2.rightJustify(leng);
        }
        
        return str2;
    }
    
    
    /*将C语言的科学计数法转换为fortran语言的表示方式，主要用于输出*/
    GString GString::sci2for(  const GString::size_type startPos,
                             const GString::size_type length,
                             const GString::size_type expLen,
                             const bool checkSwitch
                             )
    throw(GStringException)
    {
        try
        {
            GString::size_type idx = find('.', startPos);
            int expAdd = 0;
            GString exp;
            long iexp;
            //If checkSwitch is false, always redo the exponential. Otherwise,
            //set it to false.
            bool redoexp = !checkSwitch;
            
            // Check for decimal place within specified boundaries
            if ((idx == 0) || (idx >= (startPos + length - expLen - 1)))
            {
                GStringException e("sci2for: no decimal point in string");
                GFC_THROW(e);
            }
            
            // Here, account for the possibility that there are
            // no numbers to the left of the decimal, but do not
            // account for the possibility of non-scientific
            // notation (more than one digit to the left of the
            // decimal)
            if (idx > startPos)
            {
                redoexp = true;
                //const_iterator it = begin() + idx;
                // Swap digit and decimal.
                *(begin() + idx ) = *(begin()+idx-1);
                *(begin()+idx -1) = '.';
                // Only add one to the exponent if the number is non-zero
                if (  substr(startPos, length).asDOUBLE() != 0.0)
                    expAdd = 1;
            }
            
            idx =  find('e', startPos);
            if (idx == GString::npos)
            {
                idx = find('E', startPos);
                if (idx == GString::npos)
                {
                    GStringException e("sci2for:no 'e' or 'E' in string");
                    GFC_THROW(e);
                }
            }
            // Change the exponent character to D normally, or E of checkSwitch is false.
            if ( checkSwitch )
            {
                *(begin()+idx) ='D';
                //aStr[idx] = 'D';
            }
            else
            {
                *(begin()+idx) ='E';
                //aStr[idx] = 'E';
            }
            // Change the exponent itself
            if (redoexp)
            {
                exp = substr(idx + 1, GString::npos);
                iexp = exp.asINT();
                iexp += expAdd;
                
                erase(idx + 1);
                if (iexp < 0)
                {
                    *this += "-";
                    iexp -= iexp*2;
                }
                else
                {
                    *this += "+";
                }
                GString tmp(iexp) ;
                *this +=   tmp.rightJustify(expLen,'0');
            }
            
            // if the number is positive, append a space
            // (if it's negative, there's a leading '-'
            if( *(begin()) == '.')
            {
                insert((GString::size_type)0, 1, ' ');
            }
            
            //If checkSwitch is false, add on one leading zero to the string
            if ( !checkSwitch )
            {
                insert((GString::size_type)1, 1, '0');
            }
            
            return *this;
        }
        catch(GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }  // end sci2for
    
    
    
    GString GString::doub2for( const GString::size_type length,
                              const GString::size_type expLen,
                              const bool checkSwitch)
    throw(GStringException)
    {
        try
        {
            short exponentLength = expLen;
            
            /* Validate the assumptions regarding the input arguments */
            if (exponentLength < 0) exponentLength = 1;
            if (exponentLength > 3 && checkSwitch) exponentLength = 3;
            
            GString str1 = double2Sci(length, exponentLength, true, checkSwitch);
            GString toReturn = str1.sci2for( 0, length, exponentLength, checkSwitch);
            
            return toReturn;
        }
        catch( GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch( std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    
    double GString::for2doub( const GString::size_type startPos,
                             const GString::size_type length)
    {
        GString s(*this, startPos, length);  //首先 截取子串
        strip_v(s);  //剔除空格
        
        // you can blame Rinex for these special checks
        if (s.empty())
        {
            return 0;
        }
        
        GString::size_type pos = s.find_first_of("EDd");
        if ( pos != GString::npos)
        {
            s[pos] = 'e';
        }
        else
        {
            // just treat it like a double
            return substr(startPos,length).asDOUBLE();
            //return asDouble(aStr.substr(startPos, length));
        }
        
        std::stringstream st;
        st << s;
        
        double d;
        st >> d;
        
        return d;
    }
    
    
    GString GString::printable() throw(GStringException)
    {
        try
        {
            GString rv(*this);
            
            for ( int i = 0; i < (int)rv.length(); i++)
            {
                char c = rv[i];
                if ( !isprint(c))
                {
                    if ( iscntrl(c))
                    {
                        rv.replace(i,1,2,'^');
                        rv.replace(i+1,1,1, 64+(c));
                    }
                    else
                    {
                        GString mess( rv.substr(i,1).c2x() );
                        //GString mess( c2x(rv.substr(i,1)) );
                        rv.replace(i,1,"<"+mess+">");
                    }
                }
            }
            
            
            return rv;
        }
        catch(GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + std::string(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    
    GString GString::prettyPrint( const GString& lineDelim,
                                 const GString& indent,
                                 const GString& firstIndent,
                                 const GString::size_type len,
                                 const char wordDelim)
    throw(GStringException)
    {
        try
        {
            // chop aStr up into words based on wordDelim
            std::list<GString> wordList;
            GString tempStr(*this);
            tempStr.stripLeading_v(wordDelim);
            while (!tempStr.empty())
            {
                GString theFirstWord = tempStr.firstWord(wordDelim);
                wordList.push_back(theFirstWord);
                tempStr.stripLeading_v(theFirstWord);
                tempStr.stripLeading_v( wordDelim);
            }
            
            // now reassemble the words into sentences
            GString toReturn;
            GString thisLine = firstIndent, lastLine;
            while (!wordList.empty())
            {
                lastLine = thisLine;
                if ( !lastLine.empty())
                {
                    thisLine += wordDelim;
                }
                
                thisLine += wordList.front();
                
                if (thisLine.length() > len)
                {
                    // if the first word is longer than a line, just add it.
                    // if this is the first line, remember to add the indent.
                    if (lastLine.empty())
                    {
                        if (toReturn.empty())
                        {
                            lastLine += firstIndent;
                        }
                        
                        lastLine = wordList.front();
                    }
                    
                    toReturn += lastLine + lineDelim;
                    
                    thisLine.erase();
                    lastLine.erase();
                    
                    thisLine = indent;
                }
                else
                {
                    wordList.erase(wordList.begin());
                }
            }
            
            if ( !thisLine.empty())
            {
                toReturn += (thisLine + lineDelim);
            }
            
            *this = toReturn;
            return *this;
        }
        catch(GStringException &e)
        {
            GFC_RETHROW(e);
        }
        catch(std::exception &e)
        {
            GStringException strexc("Exception thrown: " + GString(e.what()));
            GFC_THROW(strexc);
        }
    }
    
    
}
