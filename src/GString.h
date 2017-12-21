//
//  GString.h
//  GFC
//
//  Created by lizhen on 15/9/14.
//  Copyright (c) 2015年 lizhen. All rights reserved.
//


#ifndef GFC_GString_h
#define GFC_GString_h

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


#include "Platform.h"
#include "GException.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <vector>
//#include <cstdio>   /// @todo Get rid of the stdio.h dependency if possible.
#include <cctype>
#include <limits>

#ifdef _WIN32
#if _MSC_VER < 1700
// For lower version of visual studio 2012 use gnu regex(版本为2.7)
#include "../gnuRegexLib/regex.h"
#pragma comment(lib, "../gnuRegexLib/regex.lib")
#else
// visual studio 2012 support c++ 0x, and we use std::regex
#include <regex>
#endif
#else
// TODO: we should use std::regex for upper than g++ 4.6
#include <regex.h>
#endif

namespace gfc
{
    
    /// This is thrown instread of a std::exception when a
    /// GFC::GString function fails.
    /// @ingroup exceptiongroup
    NEW_GEXCEPTION_CLASS(GStringException, GException);
    
    /// Class for configuring the appearance of hexDumpData() output
    class HexDumpDataConfig
    {
    public:
        HexDumpDataConfig()
        : showIndex(true), hexIndex(true), upperHex(false),
        idxDigits(4), indexWS(1), groupBy(1), groupWS(1),
        group2By(8), group2WS(2), bytesPerLine(16), showText(true),
        separator(0), textWS(4)
        {}
        HexDumpDataConfig(bool ashowIndex, bool ahexIndex, bool aupperHex,
                          unsigned aidxDigits, unsigned aindexWS,
                          unsigned agroupBy, unsigned agroupWS,
                          unsigned agroup2By, unsigned agroup2WS,
                          unsigned abytesPerLine, bool ashowText,
                          char aseparator, unsigned atextWS)
        : showIndex(ashowIndex), hexIndex(ahexIndex),
        upperHex(aupperHex), idxDigits(aidxDigits),
        indexWS(aindexWS), groupBy(agroupBy), groupWS(agroupWS),
        group2By(agroup2By), group2WS(agroup2WS),
        bytesPerLine(abytesPerLine), showText(ashowText),
        separator(aseparator), textWS(atextWS)
        {}
        bool showIndex; ///< display index into string on each line.
        bool hexIndex; ///< if true, use hex index numbers (else decimal).
        bool upperHex; ///< if true, use upper-case hex digits.
        unsigned idxDigits; ///< number of positions to use for index.
        unsigned indexWS; ///< number of whitespace charaters between index and data.
        unsigned groupBy; ///< number of bytes of data to show between spaces.
        unsigned groupWS; ///< number of whitespace charaters between groups of hex data.
        unsigned group2By; ///< number of groups to show per 2nd layer group (0=none, must be multiple of groupBy).
        unsigned group2WS; ///< number of whitespace charaters between 2nd layer groups.
        unsigned bytesPerLine; ///< number of bytes to display on a line of output (must be evenly divisible by both groupBy and group2By).
        bool showText; ///< if true, show text of message (unprintable characters become '.'.
        char separator; ///< character to offset text with (0 = none).
        unsigned textWS; ///< number of whitespace characters between hex and text.
    };

    
    
    //@类的定义及实现从string类继承而来，因此包含string全部的功能，是string类的一个超集
    class GString : public std::string
    {
        
    public:
        /**
         * Default constructor
         *
         * Constructs an empty GString ("")
         */
        GString() : std::string() { }
        
        /**
         * Duplicate the STL string copy constructor
         *
         * @param[in] s   The string to copy
         * @param[in] pos The starting position in the string to copy from
         * @param[in] n   The number of characters to copy
         */
        GString(const GString &s, size_type pos = 0, size_type n = npos) : std::string(s, pos, n) { }
        
        /**
         * Construct an GString from a null-terminated character array
         *
         * @param[in] s The character array to copy into the new string
         */
        GString(const value_type *s) : std::string(s) { }
        
        /**
         * Construct an GString from a character array and a length
         *
         * @param[in] s The character array to copy into the new string
         * @param[in] n The number of characters to copy
         */
        GString(const value_type *s, size_type n) : std::string(s, n) { }
        
        
        
        /**
         * Create an GString with @p n copies of @p c
         *
         * @param[in] n The number of copies
         * @param[in] c The character to copy @p n times
         */
        GString(size_type n, value_type c) : std::string(n, c) { }
        
        
        /*
          *  新增加的字符串构造函数
          *  由long double数据直接构造字符串
          */
        GString( const long double x, const GString::size_type precision = 21)
        {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision( static_cast<int>(precision) ) << x ;
            std::string mys = ss.str();
            GString mygs(mys.data());
            *this = mygs;
        }
        
        /*
         *  新增加的字符串构造函数
         *  由 double数据直接构造字符串
         */
        GString( const double x, const GString::size_type precision = 21)
        {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(static_cast<int_t>(precision)) << x;
            std::string mys = ss.str();
            GString mygs(mys.data());
            *this = mygs;
        }
        
        /*
         *  新增加的字符串构造函数
         *  由 int数据直接构造字符串
         */
        GString( int x )
        {
            std::ostringstream ss;
            //setw(3) << setfill('0');
            //ss << std::setw(width)<<std::setfill('0')<<x;
            
            ss << std::fixed << x;
            
            std::string mys = ss.str();
            GString mygs(mys.data());
            *this = mygs;
        }
        
        /*
         *  新增加的字符串构造函数
         *  由 任何已定义的类的对象直接构造字符串
         *  前提是该类X 已经重载运算符<<
         *  由此基础，任何自定义类均需要重载流操作运算符
         * Convert a value in a string to a type specified by the template
         * class.  The template class type must have stream operators
         * defined.
         * @param x object to turn into the templatized type.
         * @return the template object of \a x.
         */
        
        
        template <class X>
        GString(  X x )
        {
            std::ostringstream ss;
            ss << x;
            
            //*this = ss.str();
            //ss << std::fixed << std::setprecision( static_cast<int>(precision) ) << x ;
            std::string mys = ss.str();
            GString mygs(mys.data());
            *this = mygs;
        }
        
        
        
        /*
          * 需要重写substr
          * substr(size_type pos, sizetype n)
          *这里并不改变原始字符串this的值
         */
        GString substr( size_type pos = 0, size_type n = GString::npos)
        {
            std::string  s = *this;
            std::string ss = s.substr(pos,n);
            GString  rs(ss.c_str());
            return rs;
        }
        /*
          * 析构函数
          *
         */
        ~GString() {}
        
        /**
         * Perform a formatted hex-dump of the (potentially) binary
         * data to the given stream.
         * @param s stream to dump data to.
         * @param data data to hex-dump.
         * @param indent indents the string by that many spaces.
         * @param cfg formatting configuration.
         */
         void hexDumpData(std::ostream& s,
                                const GString& data,
                                unsigned indent = 0,
                                HexDumpDataConfig cfg = HexDumpDataConfig());
        
        /**
         * Perform a formatted hex-dump of the (potentially) binary
         * data to the given stream.
         * @param s stream to dump data to.
         * @param data data to hex-dump.
         * @param tag string to put at the beginning of each line of output.
         * @param cfg formatting configuration.
         */
         void hexDumpData(std::ostream& s,
                                const GString& data,
                                const GString& tag,
                                HexDumpDataConfig cfg = HexDumpDataConfig());
        
        /**
         * Remove a string from the beginning of another string const version.
         * Occurrences of the string \a aString appearing
         * at the beginning of the string \a s are removed.
         * @param s string to be stripped (modified).
         * @param aString string to remove.
         * @param num maximum number of occurrences to remove.
         * @throws GGStringException if there's a std::exception thrown.
         * @return a reference to \a s.
         */
        
        GString stripLeading(const GString& aString, GString::size_type num = GString::npos) throw(GStringException);
        GString stripLeading_v(const GString& aString,GString::size_type num = GString::npos) throw(GStringException);
        
        /**
         * Remove a string from the beginning of another string.
         * Occurrences of the string \a pString appearing
         * at the beginning of the string \a s are removed.
         * @param s string to be stripped (modified).
         * @param pString string to remove.
         * @param num maximum number of occurrences to remove.
         * @throws GStringException if there's a std::exception thrown.
         * @return a reference to \a s.
         */
        GString stripLeading( const char* pString,GString::size_type num = GString::npos) throw(GStringException)
        { return stripLeading( GString(pString), num); }
        
        GString stripLeading_v( const char* pString,GString::size_type num = GString::npos) throw(GStringException)
        { return stripLeading_v( GString(pString), num); }
        
        
        /**
         * Strip character(s) from the beginning of a string.
         * Occurrences of the character \a aCharacter appearing
         * at the beginning of the string \a s are removed.
         * @param s string to be stripped (modified).
         * @param aCharacter character to remove.
         * @param num maximum number of occurrences to remove.
         * @throws GStringException if there's a std::exception thrown.
         * @return a reference to \a s.
         */
        GString stripLeading( const char aCharacter,GString::size_type num = GString::npos) throw(GStringException)
        { return stripLeading(GString(1,aCharacter), num); }
        
        GString stripLeading_v( const char aCharacter,GString::size_type num = GString::npos) throw(GStringException)
        { return stripLeading_v(GString(1,aCharacter), num); }
        
        /**
         * Strip blanks from the beginning of a string.
         * Occurrences of the space character appearing
         * at the beginning of the string \a s are removed.
         * @param s string to be stripped (modified).
         * @param num maximum number of occurrences to remove.
         * @throws GStringException if there's a std::exception thrown.
         * @return a reference to \a s.
         */
        GString stripLeading( GString::size_type num = GString::npos) throw(GStringException)
        { return stripLeading(GString(1,' '),num); }
        
        GString stripLeading_v( GString::size_type num = GString::npos) throw(GStringException)
        { return stripLeading_v(GString(1,' '),num); }
        
        
        /**
         * Remove a string from the end of another string.
         * Occurrences of the string \a aString appearing
         * at the end of the string \a s are removed.
         * @param s string to be stripped (modified).
         * @param aString string to remove.
         * @param num maximum number of occurrences to remove.
         * @throws GStringException if there's a std::exception thrown.
         * @return a reference to \a s.
         */
        
        GString stripTrailing( const GString& aString, GString::size_type num = GString::npos) throw(GStringException);
        GString stripTrailing_v( const GString& aString, GString::size_type num = GString::npos) throw(GStringException);
        
        /**
         * Remove a string from the end of another string.
         * Occurrences of the string \a pString appearing
         * at the end of the string \a s are removed.
         * @param s string to be stripped (modified).
         * @param pString string to remove.
         * @param num maximum number of occurrences to remove.
         * @throws GStringException if there's a std::exception thrown.
         * @return a reference to \a s.
         */
        GString stripTrailing(  const char* pString, GString::size_type num = GString::npos ) throw(GStringException)
        { return stripTrailing( GString(pString), num); }
        GString stripTrailing_v(  const char* pString, GString::size_type num = GString::npos ) throw(GStringException)
        { return stripTrailing_v( GString(pString), num); }
        
        
        /**
         * Strip character(s) from the end of a string.
         * Occurrences of the character \a aCharacter appearing
         * at the end of the string \a s are removed.
         * @param s string to be stripped (modified).
         * @param aCharacter character to remove.
         * @param num maximum number of occurrences to remove.
         * @throws GStringException if there's a std::exception thrown.
         * @return a reference to \a s.
         */
        GString stripTrailing(  const char aCharacter, GString::size_type num = GString::npos ) throw(GStringException)
        { return stripTrailing( GString(1,aCharacter), num); }
        
        GString stripTrailing_v(  const char aCharacter, GString::size_type num = GString::npos ) throw(GStringException)
        { return stripTrailing_v( GString(1,aCharacter), num); }
        
        
        /**
         * Strip blanks from the end of a string.
         * Occurrences of the space character appearing
         * at the end of the string \a s are removed.
         * @param s string to be stripped (modified).
         * @param num maximum number of occurrences to remove.
         * @throws GStringException if there's a std::exception thrown.
         * @return a reference to \a s.
         */
        GString stripTrailing( GString::size_type num = GString::npos) throw(GStringException)
        { return stripTrailing( GString(1,' '),num); }
        
        GString stripTrailing_v( GString::size_type num = GString::npos) throw(GStringException)
        { return stripTrailing_v( GString(1,' '),num); }
        
        /* Remove a string from the beginning and end of another string.
        * Occurrences of the string \a aString appearing
        * at the beginning and end of the string \a s are removed.
        * @param s string to be stripped (modified).
        * @param aString string to remove.
        * @param num maximum number of occurrences to remove.
        * @throws GStringException if there's a std::exception thrown.
          * @return a reference to \a s.
          */
        GString strip( const GString& aString, GString::size_type num = GString::npos) throw(GStringException);
        GString strip_v( const GString& aString, GString::size_type num = GString::npos) throw(GStringException);
        
        /**
         * Remove a string from the beginning and end of another string.
         * Occurrences of the string \a pString appearing
         * at the beginning and end of the string \a s are removed.
         * @param s string to be stripped (modified).
         * @param pString string to remove.
         * @param num maximum number of occurrences to remove.
         * @throws GStringException if there's a std::exception thrown.
         * @return a reference to \a s.
         */
        GString strip( const char* pString, GString::size_type num = GString::npos ) throw(GStringException)
        { return strip( GString(pString), num); }
        GString strip_v( const char* pString, GString::size_type num = GString::npos ) throw(GStringException)
        { return strip_v( GString(pString), num); }
        
        /**
         * Strip character(s) from the beginning and end of a string.
         * Occurrences of the character \a aCharacter appearing
         * at the beginning and end of the string \a s are removed.
         * @param s string to be stripped (modified).
         * @param aCharacter character to remove.
         * @param num maximum number of occurrences to remove.
         * @throws GStringException if there's a std::exception thrown.
         * @return a reference to \a s.
         */
        GString strip(const char aCharacter, GString::size_type num = GString::npos ) throw(GStringException)
        { return strip( GString(1,aCharacter), num); }
        
        GString strip_v(const char aCharacter, GString::size_type num = GString::npos ) throw(GStringException)
        { return strip_v( GString(1,aCharacter), num); }
        
        /* Strip blanks from the beginning and end of a string.
        * Occurrences of the space character appearing
        * at the beginning and end of the string \a s are removed.
        * @param s string to be stripped (modified).
        * @param num maximum number of occurrences to remove.
        * @throws GStringException if there's a std::exception thrown.
          * @return a reference to \a s.
          */
        GString strip(GString& s, GString::size_type num = GString::npos) throw(GStringException)
        { return strip( GString(1, ' '), num); }
        GString strip_v(GString& s, GString::size_type num = GString::npos) throw(GStringException)
        { return strip_v( GString(1, ' '), num); }
        
        /**
         * Converts all of the receiver's characters that are in the
         * first specified string to the corresponding character in
         * the second specified string.
         * @param aString string to perform translation on.
         * @param inputChars characters in \a aString to translate from.
         * @param outputChars characters to translate to.
         * @param pad pad character in the event inputChars and
         * outputChars are not equal length.  The pad character will
         * become the translated character.
         */
        GString translate( const GString& inputChars,
                          const GString& outputChars,
                          const char pad);
        
        /**
         * Changes occurrences of a specified pattern to a specified
         * replacement string.  You can specify the number of changes
         * to perform.  The default is to change all occurrences of
         * the pattern. You can also specify the position in the
         * receiver at which to begin.
         * @param aString string to perform translation on.
         * @param inputString The pattern string as a reference to an
         *   object of type string.  The library searches for the
         *   pattern string within the receiver's data.
         * @param outputString The replacement string as a reference
         *   to an object of type string. It replaces the occurrences
         *   of the pattern string in the receiver's data.
         * @param startPos The position to start the search at within
         *   the receiver's data.  The default is 0.
         * @param numChanges the number of patterns to search for and
         *   change.  The default is to change all occurrences of the
         *   pattern.
         * 本质上就是将字符串中的inputString替换为outputstring
         */
        GString Greplace( const GString&  originString,
                         const GString&   replaceString,
                         GString::size_type startPos = 0,
                         unsigned numChanges = (std::numeric_limits<unsigned>::max)() );
        
        /**
         * Right-justifies the receiver in a string of the specified
         * length. If the receiver's data is shorter than the
         * requested length (\a length), it is padded on the left with
         * the pad character (\a pad). The default pad
         * character is a blank.
         * @param s string to be modified.
         * @param length new desired length of string.
         * @param pad character to pad string with (blank by default).
         * @throws GStringException if there's a std::exception thrown.
         * @return a reference to \a s.  
         *
         * 即从右端开始截取字符串
         */
        GString rightJustify(   const GString::size_type length, const char pad = ' ')  throw(GStringException);
        GString leftJustify(    const GString::size_type length,const char pad = ' ')   throw(GStringException);
        
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
         * @throws GStringException if there's a std::exception thrown.
         * @return a reference to \a s.
         */
        GString center( const GString::size_type length, const char pad = ' ' ) throw(GStringException);
        
        /*将当前GString转换为double数据*/
        double asDOUBLE()
        { return strtod( c_str(), 0); }
        
        /**
         * Convert a string to an integer.
         * @param s string containing a number.
         * @return long integer representation of string.
         */
        long asINT()
        { return strtol( c_str(), 0, 10); }
        
        /**
         * Convert a string to an unsigned integer.
         * @param s string containing a number.
         * @return unsigned long integer representation of string.
         */
        unsigned long asUINT()
        { return strtoul( c_str(), 0, 10); }
        
        float asFLOAT() throw(GStringException);
        
        long double asLDOUBLE() throw(GStringException);
        
        /*将十进制的字符串转换为十六进制的字符串*/
        GString d2x() throw(GStringException);
        /*将十六进制字符串转换为十进制字符串*/
        GString x2d()  throw(GStringException);
        
        GString c2x() throw(GStringException);
        
        
        /*将字符串转换为大写字母*/
        GString   upperCase();
        /*将字符串变为小写字母*/
        GString  lowerCase();
        
        bool isDigitString();
        bool isDecimalString();
        bool isScientificString();
        bool isAlphaString();
        /*
          * 正则表达式匹配
         * aPattern  匹配模板 This is a POSIX regular expression.
         * zeroOrMore(*)  oneOrMore(+) anyChar(.)处理其他的情况
         *
         */
        GString matches( const GString& aPattern, const char zeroOrMore='*', const char oneOrMore='+', const char anyChar='.') throw(GStringException);
        
        /**
         * Perform pattern matching on strings.
         * Looks for a pattern in a string.  Wildcards are allowed.
         * Uses POSIX regular expressions.
         * @param s string to search.
         * @param aPattern pattern to search for. This is a POSIX
         * regular expression.
         * @param zeroOrMore character representing wildcards
         * matching strings of zero or more characters (default '*').
         * @param oneOrMore character representing plus sign
         * matching strings of one or more characters (default '+').
         * @param anyChar character representing wildcards matching a
         * single arbitrary character (default '.').
         * @return t if a match is found, f if not.
         */
        bool isLike( const GString& aPattern,const char zeroOrMore = '*',const char oneOrMore = '+',const char anyChar = '.' )
        throw(GStringException)
        { return matches(aPattern, zeroOrMore, oneOrMore, anyChar) != GString(); }
        
        
        /**
         * Work-horse method for printf.  Substitutes patterns
         * matching \a pat with \a rep.  Use only one pattern/token
         * at a time!  This used to be DayTime::iprint().
         * @param fmt format to use for this time.
         * @param pat regular expression pattern to match.
         * @param rep sprintf token replacement.  First character is
         * token character used in fmt, remainder is sprintf token to
         * use.  For example, with fmt="%15S", pat="%[ 0-]?[[:digit:]]*S",
         * and rep="Sd", the fmt will be translated to "%15d" before
         * using it in a sprintf call like printf("%15d"), \a to.
         * @param to the value to stuff into the string.
         * @return \a fmt with \a pat replaced by \a to.  If there is no
         * match, \a fmt is returned unchanged.
         */
        template <class T>
        GString formattedPrint( GString& fmt,  GString& pat,  GString& rep, T to ) throw(GStringException);
        
       /**
         * Returns the first word in string \a s without modifying the string.
         * @param s the string to count the words from.
         * @param delimiter the character that marks the start and
         * end of a word.
         * @return the first word from \a s;
         */
       GString firstWord( const char delimiter = ' ' )  throw(GStringException);
        
       /*统计该字符串中有多少个word*/
       int  wordsNum(const char delimiter = ' ') throw(GStringException);
        /*返回一个子串，从第firsWord往后numWords*/
       GString wordsString( const GString::size_type firstWord = 0,const GString::size_type numWords= GString::npos,const char delimiter = ' ')
       throw(GStringException);
       
        /**
         * 取得第wordNum个word的字符串
         * Returns word number \a wordNum from \a s (if any).
         * @param s a string with the word you want removed.
         * @param wordNum the number of the word you want from \a s.
         * The first word is word 0.
         * @param delimiter the character that marks the start and
         * end of a word.
         * @return the first word from \a s or an empty string if there is
         * no \a wordNum'th word.
         */
        GString word( const GString::size_type wordNum = 0, const char delimiter = ' ')
        throw(GStringException)
        { return wordsString(wordNum, 1, delimiter); }
        
        
        /*去掉第一个word；返回去掉后的子串*/
        GString stripFirstWord(const char delimiter = ' ') throw(GStringException);
        GString stripFirstWord_v(const char delimiter = ' ') throw(GStringException);
        
        /*字符串分割*/
        std::vector<GString> split(const char delimiter = ' ') throw(GStringException);
        
        /*
          *删除从first开始往后wordsToReplace个数的words
         * 改变本身的值
         * Removes indicated words from the string \a s.
         * \a s is modified as a result.
         * @param s a string with the words you want removed.
         * @param first the first word to be removed (the first word is 0).
         * @param wordsToReplace the number of words you want removed.
         * @param delimiter the character that marks the start and
         * end of a word.
         * @return a reference to string \a s with the words removed.
         */
       GString removeWords( const GString::size_type first = 0 ,const GString::size_type wordsToReplace = GString::npos,const char delimiter= ' ')
        throw(GStringException);
        
       
        /** Convert a double GString to scientific notation; this routine works better,
         * on Windows particularly, than doub2sci.
         * @param length = total string length,
         *                         including 1 for overall sign if showPlus is true.
         * @param precision = number of digits after the decimal and before the 'e'
         * @param explen = length of exponent, this must = 1, 2 or 3; 最大支持999次方
         * NB. length is increased if precision, explen and showPlus require it.
         */
        GString double2Sci( const GString::size_type length,
                                     const GString::size_type precision,
                                     const GString::size_type explen,
                                     bool showPlus = false);
        
        
        /**
         * Convert scientific notation to FORTRAN notation.
         * As an example, the string "1.5636E5" becomes " .15636D6".
         * Note that the first character of the string will be '-' if
         * the number is negative or ' ' if the first character is positive.
         * @param aStr string with number to convert
         * @param startPos start position of number in string
         * @param length length (in characters) of number, including exponent.
         * @param expLen length (in characters of exponent, not including sign.
         * @param checkSwitch will keep the method running as orignially programed
         * when set to true.  If false, the method will always resize exponentials,
         * produce an exponential with an E instead of a D, and always have a leading
         * zero.  For example -> 0.87654E-0004 or -0.1234E00005.
         * @throws Exception if the string is not a number in scientific notation
         */
        
        GString sci2for( const GString::size_type startPos = 0,
                        const GString::size_type length = GString::npos,
                        const GString::size_type expLen = 3,
                        const bool checkSwitch = true
                                 )
        throw(GStringException);
        
        /**
         * Convert double precision floating point to a string
         * containing the number in FORTRAN notation.
         * As an example, the number 156360 becomes ".15636D6".
         * @param d number to convert.
         * @param length length (in characters) of number, including exponent.
         * @param expLen length (in characters of exponent, including sign.
         * @param checkSwitch if true, keeps the exponential sanity check for
         * exponentials above three characters in length.  If false, it removes
         * that check.
         * @return a string containing \a d in FORTRAN notation.
         */
        GString doub2for( const GString::size_type length,
                        const GString::size_type expLen,
                        const bool checkSwitch)
        throw(GStringException);
        
        
        /**
         * Convert FORTRAN representation of a double precision
         * floating point in a string to a number.
         * As an example, the number ".15636D6" becomes 156360.
         * @param aStr string containing FORTRAN representation of number.
         * @param startPos beginning of number in string.
         * @param length length (in characters) of number, including exponent.
         * @return value of the number.
         * 用于读取fortran程序输出的数据文档
         */
        double for2doub( const GString::size_type startPos = 0, const GString::size_type length = GString::npos);
        
        
        /**
         * Change a string into printable characters.  Control
         * characters (0-26) are changed to ^@, ^A, etc.  Other
         * non-printable characters are changed to hex sequences
         * enclosed in <>.
         * @param aStr the string to make printable.
         */
       GString printable() throw(GStringException);
        
       
        /**
         *  专门为RINEX文件读写打造的，哈哈哈哈
         * Const version of prettyPrint, which nicely expands the
         * input string into several lines.
         * @param aStr the string to be modified.
         * @param lineDelim a string to put between every line.
         * @param indent an indentataion string used on all but the first line
         * @param firstIndent is the indentation used on the first line.
         * @param len the maximum length of string to put on a line.
         * @param wordDelim the character that separates each word.
         * @return the string nicely formatted.
         */
      GString prettyPrint( const GString& lineDelim ="\n",
                                     const GString& indent= "",
                                     const GString& firstIndent = "     ",
                                     const GString::size_type len = 80 ,
                                     const char wordDelim = ' ')
        throw(GStringException);
        
        
        /**
         * Convert a value in a string to a type specified by the template
         * class.  The template class type must have stream operators
         * defined.
         * @param x object to turn into the templatized type.
         * @return the template object of \a x.
         */
        template <class X>
        X asData() throw(GStringException)
        {
            try
            {
                GString s = *this;
                std::istringstream is(s);
                X x;
                is >> x;
                return x;
            }
            catch(std::exception &e)
            {
                GStringException strexc("Exception thrown: " + GString(e.what()));
                GFC_THROW(strexc);
            }
        }
        
        
        template<class X>
         GString asString(const X x)
        {
            std::ostringstream ss;
            ss << x;
            return ss.str();
        }
        
        
        
        
    private:
        //没有私有成员变量
        
    };  // end of the class GString
    
    
    
    
    
    //成员函数具体实现
    
}  // end of namespace

#endif  //end of GFC_GString_H
