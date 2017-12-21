
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


/****************************************************************************
purpose：   used for the process of JPL planet ephemeris
time：2014.11.08
author：    Zhen.Li
version:     V1.0
    All Rights Reserved
****************************************************************************/

#include "GJPLEPH.h"

namespace gfc {

// the definition of static variables
GJPLEPH::JPLHEADER GJPLEPH::ephHeader;
std::vector<GJPLEPH::JPLRECORD> GJPLEPH::ephData;
double GJPLEPH::m_AU = 0.0;  // astronomial unit,
double GJPLEPH::m_EMRAT = 0.0;

GJPLEPH::GJPLEPH() {
  m_bary = false;
  m_km = false;

  // int h = static_cast<int>(ephHeader.ephindicator.size() );
  // if( h == 0 )
  //{
  //    printf("ERROR: please reading the planet eph data file first!\n");
  //}
}

GJPLEPH::~GJPLEPH() {}

/*KM = 1,单位为km，km = 0 单位为AU*/
void GJPLEPH::SetUnit(bool km) { m_km = km; }

/*
BARY = 0 为太阳系质心, BARY = 1为 太阳质心
*/
void GJPLEPH::SetBARY(bool bary) { m_bary = bary; }

void GJPLEPH::loadEphFile_a(GString filename) {
  // open the file
  ifstream strm;
  strm.open(filename);
  if (!strm) {
    GException e("Failed to open input file " + filename + ". Abort.");
    GFC_THROW(e);
  }

  GString line, word;
  int group = 0, n = 0;  // n will count lines/items within a group
  int m = 0;
  std::vector<GString> const_names;
  while (1) {
    getline(strm, line);
    line.stripTrailing_v("/r");

    // catch new groups
    if (line.substr(0, 5) == "GROUP") {
      word = line.stripFirstWord();  // stripFirstWord(line);
      group = static_cast<int>(
          line.stripFirstWord().asINT());  //  stripFirstWord(line);
      n = 0;  // n will count lines/items within a group
      m = 0;
      continue;
    }

    // skip blank lines
    line.stripLeading_v(" ");
    if (line.empty()) {
      if (strm.eof() || !strm.good())
        break;  // if the last line is blank
      else
        continue;
    }

    // process entire line at once
    // first line (no GROUP)
    if (group == 0) {
      word = line.stripFirstWord_v();  // stripFirstWord(line);
      if (word == "KSIZE=") {
        ephHeader.KSIZE = static_cast<int>(line.stripFirstWord_v().asINT());
      }
      word = line.stripFirstWord_v();  // stripFirstWord(line);
      if (word == "NCOEFF=") {
        ephHeader.NCOEFF = static_cast<int>(
            line.firstWord().asINT());  // asInt(stripFirstWord(line));
        continue;
      } else {
        GException e("Confused on the first line - 3rd word is not NCOEFF=");
        GFC_THROW(e);
      }
    }
    // GROUP 1010
    else if (group == 1010) {
      if (n > 2) {  // this should not happen
        GException e("Too many labels under GROUP 1010");
        GFC_THROW(e);
      } else {
        n++;
        continue;
      }
    }
    // GROUP 1030
    else if (group == 1030) {
      // start and stop times. These are meaningless here, because they will be
      // determined by the data that follows this header, and so are meaningful
      // only in the binary file.
      ephHeader.start_jed = line.stripFirstWord_v().for2doub();
      ephHeader.end_jed = line.stripFirstWord_v().for2doub();
      ephHeader.record_span = line.stripFirstWord_v().for2doub();
    }
    // GROUP 1070 - end-of-header and the start of the record
    else if (group == 1070) {
      int iret = 0;
      // expect this many lines per record
      int nmax = ephHeader.NCOEFF / 3 + (ephHeader.NCOEFF % 3 ? 1 : 0);
      // loop over lines in the file
      int ntot = 0;  // counts the total number of lines
      int n = 0;     // counts the lines within a set of coefficients
      int nc = 0;    // count coefficients within a record
      std::vector<double> data_vector;
      JPLRECORD testrecord;
      n = 1;
      testrecord.index = static_cast<int>(line.stripFirstWord_v().asINT());
      testrecord.data_num = static_cast<int>(line.stripFirstWord_v().asINT());

      while (1) {
        getline(strm, line);
        line.stripTrailing_v("\r");
        if (line.empty()) {
          if (strm.eof()) break;
          if (!strm.good()) {
            iret = -1;
            break;
          }
          continue;
        }

        if (n == 0)  // start a new record
        {
          testrecord.index = static_cast<int>(line.stripFirstWord_v().asINT());
          int ncc =
              static_cast<int>(line.stripFirstWord_v()
                                   .asINT());  //  asInt(stripFirstWord(line));
                                               //  // 2nd word is ncoeff
          testrecord.data_num = ncc;

          if (ncc != ephHeader.NCOEFF) {
            GException e("readASCIIdata finds conflicting sizes in header (" +
                         GString(ephHeader.NCOEFF) + ") and data (" +
                         GString(ncc) + ") in file " + filename + " at line #" +
                         GString(ntot));
            GFC_THROW(e);
          }
          nc = 0;
        } else {
          for (int j = 0; j < 3; j++) {
            double coeff = line.stripFirstWord_v()
                               .for2doub();  // for2doub(stripFirstWord(line));
            nc++;
            testrecord.record.push_back(coeff);
            if (nc >= ephHeader.NCOEFF)  // the end of the record
            {
              ephData.push_back(testrecord);
              testrecord.record.clear();
              break;
            }
          }
        }

        if (strm.eof()) break;
        if (!strm.good()) {
          iret = -1;
          break;
        }
        if (n == nmax)
          n = 0;
        else
          n++;
        ntot++;
      }

      break;
    }

    // process the line one (whitespace-separated) word at a time
    while (!line.empty()) {
      word = line.stripFirstWord_v();  // stripFirstWord(line);

      if (group == 1040) {
        if (n++ == 0) {
          ephHeader.const_num = static_cast<int>(word.asINT());
        } else {
          const_names.push_back(word);
        }
      } else if (group == 1041) {
        if (n++ == 0) {
          if ((int)ephHeader.const_num != word.asINT()) {
            GString mystr;
            mystr.asString(ephHeader.const_num);
            GException e("Nconst does not match N in GROUP 1041 : " + mystr +
                         " != " + word);
            GFC_THROW(e);
          }
        } else {
          ephHeader.constant[const_names[n - 2]] = word.for2doub();
        }
      } else if (group == 1050) {
        if (m == 0) {
          INDICATOR myindicator;
          myindicator.start_pos = static_cast<int>(word.asINT());
          ephHeader.ephindicator.push_back(myindicator);
        } else if (m == 1) {
          ephHeader.ephindicator[n - ephHeader.ephindicator.size()].num_coeff =
              static_cast<int>(word.asINT());
        } else if (m == 2) {
          ephHeader.ephindicator[n - 2 * ephHeader.ephindicator.size()]
              .num_blocks = static_cast<int>(word.asINT());
        }
        n++;
      } else {
        GString mystr;
        mystr.asString(group);
        GException e("Confused about GROUP : " + mystr);
        GFC_THROW(e);
      }
    }  // end loop over words

    m++;
    if (strm.eof() || !strm.good()) break;  // if the last line is not blank
  }

  strm.clear();
  strm.close();

  m_AU = ephHeader.constant["AU"];             //天文单位常量,单位为km
  
  //GFCCONST::RegByName("AU", m_AU);
    
  m_EMRAT = ephHeader.constant["EMRAT"];       //地月质量比
  ephHeader.start_jed = ephData[0].record[0];  //第一个数据记录的开始时刻
  ephHeader.end_jed =
      ephData[ephData.size() - 1].record[1];  //最后一个数据记录的结束时刻
  ephHeader.record_span =
      ephData[0].record[1] -
      ephData[0].record[0];  //时间间隔由第一个数据记录计算出来
  //
}

/*
 *
 * all the time should be in TDB
 *
 */
void GJPLEPH::setEpoch(GTime et) {
  // m_epoch = epoch;
  // need to get et2[2] from m_epoch

  m_epoch = et;
    
  double et2[2] = {0.0}, pjd[4] = {0.0}, s = 0.0;
  JDTime jdt = GTime::GTime2JDTime(et);

  et2[0] = jdt.m_jd;
  et2[1] = (jdt.m_sod + jdt.m_fsod) / 86400.0;

  s = et2[0] - 0.5;
  datasplit(s, &pjd[0]);
  datasplit(et2[1], &pjd[2]);
  pjd[0] = pjd[0] + pjd[2] + 0.5;
  pjd[1] = pjd[1] + pjd[3];
  datasplit(pjd[1], &pjd[2]);
  pjd[0] = pjd[0] + pjd[2];
  /* here pjd[0] contains last midnight before epoch desired (in JED: *.5)
   and pjd[3] contains the remaining, fractional part of the epoch         */
  /*   error return for epoch out of range  */
  if ((pjd[0] + pjd[3]) < ephHeader.start_jed ||
      (pjd[0] + pjd[3]) > ephHeader.end_jed)
  {
    puts("Requested Epoch time is NOT in ephemeris limits.\n");
    return;
  }

  /*   calculate record # and relative time in interval   */
  m_nr = (long)((pjd[0] - ephHeader.start_jed) / ephHeader.record_span);
  /* add 2 to adjust for the first two records containing header data
   * 修改原来有+2 ，现在直接去掉*/
  if (pjd[0] == ephHeader.end_jed)
  {
    m_nr = m_nr - 1;
  }
    
  m_t[0] =
      (pjd[0] - ((1.0 * m_nr) * ephHeader.record_span + ephHeader.start_jed) +
       pjd[3]) /
      ephHeader.record_span;
  if (m_km) {
    m_t[1] = ephHeader.record_span * 86400.0;  // velocity 单位换算为m/s
    m_aufac = 1.0;
  } else {
    m_t[1] = ephHeader.record_span;  // velocity 单位换算为 AU/day
    m_aufac = 1.0 / m_AU;
  }
    
  double pefau[6] = {0.0};
    
  int index = GJPLEPH::SUN - 1;  // the 11th is sun
  int pos = ephHeader.ephindicator[index].start_pos;
  int coeff = ephHeader.ephindicator[index].num_coeff;
  int blocks = ephHeader.ephindicator[index].num_blocks;
    
  /*  every time interpolate Solar System barycentric sun state   */
  interp(&ephData[m_nr].record[pos - 1], m_t, coeff, 3, blocks, 2, pefau);

  for (int i = 0; i < 6; ++i) {
    m_sun_SB[i] = pefau[i] * m_aufac;
  }
}

void GJPLEPH::setEpoch(double et) {
  // m_epoch = epoch;
  // need to get et2[2] from m_epoch

  double et2[2] = {0.0}, pjd[4] = {0.0}, s = 0.0;
  et2[0] = et;
  et2[1] = 0.0;

  s = et2[0] - 0.5;
  datasplit(s, &pjd[0]);
  datasplit(et2[1], &pjd[2]);
  pjd[0] = pjd[0] + pjd[2] + 0.5;
  pjd[1] = pjd[1] + pjd[3];
  datasplit(pjd[1], &pjd[2]);
  pjd[0] = pjd[0] + pjd[2];
  /* here pjd[0] contains last midnight before epoch desired (in JED: *.5)
   and pjd[3] contains the remaining, fractional part of the epoch         */
  /*   error return for epoch out of range  */
  if ((pjd[0] + pjd[3]) < ephHeader.start_jed ||
      (pjd[0] + pjd[3]) > ephHeader.end_jed) {
    puts("Requested Epoch time is NOT in ephemeris limits.\n");
    return;
  }

  /*   calculate record # and relative time in interval   */
  m_nr = (long)((pjd[0] - ephHeader.start_jed) / ephHeader.record_span);
  /* add 2 to adjust for the first two records containing header data
   * 修改原来有+2 ，现在直接去掉*/
  if (pjd[0] == ephHeader.end_jed) {
    m_nr = m_nr - 1;
  }

  m_t[0] =
      (pjd[0] - ((1.0 * m_nr) * ephHeader.record_span + ephHeader.start_jed) +
       pjd[3]) /
      ephHeader.record_span;
  if (m_km) {
    m_t[1] = ephHeader.record_span * 86400.0;  // velocity 单位换算为m/s
    m_aufac = 1.0;
  } else {
    m_t[1] = ephHeader.record_span;  // velocity 单位换算为 AU/day
    m_aufac = 1.0 / m_AU;
  }

  double pefau[6] = {0.0};

  int index = GJPLEPH::SUN - 1;  // the 11th is sun
  int pos = ephHeader.ephindicator[index].start_pos;
  int coeff = ephHeader.ephindicator[index].num_coeff;
  int blocks = ephHeader.ephindicator[index].num_blocks;

  /*  every time interpolate Solar System barycentric sun state   */
  interp(&ephData[m_nr].record[pos - 1], m_t, coeff, 3, blocks, 2, pefau);

  for (int i = 0; i < 6; ++i) {
    m_sun_SB[i] = pefau[i] * m_aufac;
  }
}

/****************************************************************************/
/*****************************************************************************
 **                         pleph(et,ntar,ncent,rrd)                         **
 ******************************************************************************
 **                                                                          **
 **    This subroutine reads the jpl planetary ephemeris                     **
 **    and gives the position and velocity of the point 'ntarg'              **
 **    with respect to 'ncent'.                                              **
 **                                                                          **
 **    Calling sequence parameters:                                          **
 **                                                                          **
 **      et = (double) julian ephemeris date at which interpolation          **
 **           is wanted.                                                     **
 **                                                                          **
 **    ntarg = integer number of 'target' point.                             **
 **                                                                          **
 **    ncent = integer number of center point.                               **
 **                                                                          **
 **    The numbering convention for 'ntarg' and 'ncent' is:                  **
 **                                                                          **
 **            1 = mercury           8 = neptune                             **
 **            2 = venus             9 = pluto                               **
 **            3 = earth            10 = moon                                **
 **            4 = mars             11 = sun                                 **
 **            5 = jupiter          12 = solar-system barycenter             **
 **            6 = saturn           13 = earth-moon barycenter               **
 **            7 = uranus           14 = nutations (longitude and obliq)     **
 **                                 15 = librations, if on eph. file         **
 **                                                                          **
 **            (If nutations are wanted, set ntarg = 14.                     **
 **             For librations, set ntarg = 15. set ncent= 0)                **
 **                                                                          **
 **     rrd = output 6-element, double array of position and velocity        **
 **           of point 'ntarg' relative to 'ncent'. The units are au and     **
 **           au/day. For librations the units are radians and radians       **
 **           per day. In the case of nutations the first four words of      **
 **           rrd will be set to nutations and rates, having units of        **
 **           radians and radians/day.                                       **
 **                                                                          **
 **           The option is available to have the units in km and km/sec.    **
 **           for this, set km=TRUE at the beginning of the program.         **
 *****************************************************************************/
/*
*
THE NUMBERING CONVENTION FOR 'NTARG' AND 'NCENT' IS:
C
C                1 = MERCURY           8 = NEPTUNE
C                2 = VENUS             9 = PLUTO
C                3 = EARTH            10 = MOON
C                4 = MARS             11 = SUN
C                5 = JUPITER          12 = SOLAR-SYSTEM BARYCENTER
C                6 = SATURN           13 = EARTH-MOON BARYCENTER
C                7 = URANUS           14 = NUTATIONS (LONGITUDE AND OBLIQ)
C                          15 = LIBRATIONS, IF ON EPH FILE
C
C             (IF NUTATIONS ARE WANTED, SET NTARG = 14. FOR LIBRATIONS,
 C                   SET NTARG = 15. SET NCENT=0.)
*
*
* 1         2       3                      4      5
* Mercury, Venus, Earth-Moon barycenter, Mars, Jupiter, Saturn, Uranus, Neptune,
Pluto,Moon and Sun, Nutation, Lunar libration
*
*/
void GJPLEPH::getPos(int ncent, int ntarg, double* pv) {
  memset(pv, 0, sizeof(double) * 6);

  bool bsave = false;

  if (ntarg == ncent) return;

  if (ntarg == GJPLEPH::NUTATIONS) {
    calPos(12, 2, pv);
    return;
  } else if (ntarg == GJPLEPH::LIBRATIONS) {
    calPos(13, 2, pv);
    return;
  }

  bsave = m_bary;
  m_bary = true;
    
  double pv_dummy[13][6] = {{0}};
  int list[12] = {0};
  /* list is a vector denoting, for which "body"
   ephemeris values should be calculated by state():
   0=Mercury,1=Venus,2=EMBary,...,8=Pluto,
   9=geocentric Moon, 10=nutations in long. & obliq.
   11= lunar librations  */
  int i = 0;
  /*  set up proper entries in 'list' array for state call     */
  for (i = 0; i < 2; ++i) /* list[] IS NUMBERED FROM ZERO ! */
  {
    int k = ntarg - 1;
    if (i == 1) k = ncent - 1; /* same for ntarg & ncent */
    if (k <= 9) list[k] = 2;   /* Major planets */
    if (k == 9) list[2] = 2;   /* for moon state earth state is necessary*/
    if (k == 2) list[9] = 2;   /* for earth state moon state is necessary*/
    if (k == 12) list[2] = 2;  /* EMBary state additionaly */
  }

  for (i = 0; i < 10; i++) {
    if (list[i] == 0) {
      continue;
    }
    calPos(i + 1, 2, &pv_dummy[i][0]);
  }

  /* Solar System barycentric Sun state goes to pv[10][] */
  if (ntarg == 11 || ncent == 11)
    for (i = 0; i < 6; ++i) pv_dummy[10][i] = m_sun_SB[i];

  /* Solar System Barycenter coordinates & velocities equal to zero */
  if (ntarg == 12 || ncent == 12)
    for (i = 0; i < 6; ++i) pv_dummy[11][i] = 0.0;

  /* Solar System barycentric EMBary state:  */
  if (ntarg == 13 || ncent == 13)
    for (i = 0; i < 6; ++i) pv_dummy[12][i] = pv_dummy[2][i];

  /* if moon from earth or earth from moon ..... */
  if ((ntarg * ncent) == 30 && (ntarg + ncent) == 13)
    for (i = 0; i < 6; ++i) pv_dummy[2][i] = 0.0;
  else {
    if (list[2] == 2) /* calculate earth state from EMBary */
      for (i = 0; i < 6; ++i)
        pv_dummy[2][i] -= pv_dummy[9][i] / (1.0 + m_EMRAT);

    if (list[9] == 2) /* calculate Solar System barycentric moon state */
      for (i = 0; i < 6; ++i) pv_dummy[9][i] += pv_dummy[2][i];
  }

  for (int i = 0; i < 6; i++) {
    pv[i] = pv_dummy[ntarg - 1][i] - pv_dummy[ncent - 1][i];
  }

  m_bary = bsave;
}

    
/*
*  get the apparent position,
*
*/
void GJPLEPH::getPos_apparent( int ncent, int ntarg, double* pv)
{
    GTime epoch  = m_epoch;  //should be in TDB
    double tPv[6] = {0.0};
    double c = GCONST("CLIGHT")/1000.0;
    
    double t0 = 0.0, t1 = 0.0;
    
    bool km_bake = m_km;  //use AU as unit
    m_km = false;
    
    getPos( ncent, ntarg, tPv );
    t0 = sqrt( tPv[0]*tPv[0] + tPv[1]*tPv[1] + tPv[2]*tPv[2] )/c;
    
    while (1)
    {
        SetUnit(false); // set unit AU
        
        setEpoch(epoch-t0);
        
        getPos( ncent, ntarg, tPv );
        
        //how to calculate the distance, because the light can be bent when it cross the big celescial bodies
        t1 = sqrt( tPv[0]*tPv[0] + tPv[1]*tPv[1] + tPv[2]*tPv[2] )/c*m_AU ;
        
        if( fabs(t0-t1)< 1.0E-10 )  // 0.1 ns
        {
            break;
        }
        else
        {
            t0 = t1;
        }
    }
    
    for(int i = 0 ; i< 6 ; i++ )
    {
        pv[i] = tPv[i]*m_AU;
    }
    
    
    SetUnit(km_bake);
    setEpoch(epoch);
    
    
}
    
/*
iplanet: the index of the planet, start from 1
pv:  the positiona and velocity
opt: 1 for position only, 2 for position and velocity
iplanet : start from 1 !!
order:
Mercury, Venus, Earth-Moon barycenter, Mars, Jupiter, Saturn, Uranus, Neptune,
Pluto,Moon and Sun, Nutation, Lunar libration
*/
void GJPLEPH::calPos(int iplanet, int opt, double* pv) {
  int index = iplanet - 1;
  double pefau[6] = {0.0};

  int pos = ephHeader.ephindicator[index].start_pos;
  int coeff = ephHeader.ephindicator[index].num_coeff;
  int blocks = ephHeader.ephindicator[index].num_blocks;

  if (iplanet == 12 && opt > 0 && coeff > 0)  // just for nutation
  {
    interp(&ephData[m_nr].record[pos - 1], m_t, coeff, 2, blocks, opt, pefau);
    for (int i = 0; i < 4; i++) {
      pv[i] = pefau[i];
    }
    return;
  } else if (iplanet == 13 && opt > 0 && coeff > 0)  // just for libration
  {
    interp(&ephData[m_nr].record[pos - 1], m_t, coeff, 3, blocks, opt, pefau);
    memcpy(pv, pefau, sizeof(double) * 6);
    return;
  } else {
    interp(&ephData[m_nr].record[pos - 1], m_t, coeff, 3, blocks, opt, pefau);
  }

  for (int j = 0; j < 6; ++j) {
    if (iplanet <= GJPLEPH::SUN && !m_bary)  // the first 10 (0-9) planets are
    {
      pv[j] =
          pefau[j] * m_aufac - m_sun_SB[j];  // referenced the barycenter of sun
    } else  // for the last nutaiton and libration, output the coordinated
            // referenced barycenter of solar system
    {
      pv[j] = pefau[j] * m_aufac;
    }
  }
}

void GJPLEPH::EphTest(std::string testfilename) {
  FILE* testf = NULL;
  testf = fopen(testfilename.c_str(), "r");
  if (testf == NULL) {
    printf("星历测试文件%s打开失败!\n", testfilename.c_str());
  }

  printf("line--jed--t# c# x# --jpl value---user value --- difference--\n");
  int line = 0;
  char buff[100] = {0};
  std::string strtmp;
  while (!feof(testf)) {
    fgets(buff, 100, testf);
    strtmp = buff;
    strtmp = strtmp.substr(
        0, strtmp.length() - 1);  //去掉回车符号,linux系统可能存在问题

    if (!strtmp.find("EOT"))  //跳过文件头，EOT文件头结束
    {
      double rrd[6] = {0.0};
      while (fgets(buff, 100, testf)) {
        int denum = 0, ntarget = 0, ncenter = 0, ncoord = 0;
        double jed = 0.0, xi = 0.0;
        char datetime[11] = {0};  //用于存储年月日
        sscanf(buff, "%d %s %lf %d %d %d %lf", &denum, datetime, &jed, &ntarget,
               &ncenter, &ncoord, &xi);
        if (jed > ephHeader.end_jed || jed < ephHeader.start_jed) {
          continue;
        }

        setEpoch(jed);
        getPos(ncenter, ntarget, rrd);
        // pleph(jed, ntarget, ncenter, rrd);

        double del = fabs(xi - rrd[ncoord - 1]);
        if (ntarget == 15 && ncoord == 3) {
          // JDEPOC = 2440400.5d0;
          del = del / (1.0 + 100.0 * fabs(jed - 2440400.5) / 365.25);
        }
        if (del >= 1.0E-13) {
          printf("*****  WARNING : next difference >= 1.D-13  *****\n");
          system("pause");
        }

        printf("%d %lf %d %d %d %lf %lf %e\n", line, jed, ntarget, ncenter,
               ncoord, xi, rrd[ncoord - 1], del);
        line = line + 1;
      }
    }
  }
}

/*
根据常数的名称获取常数的值
*/
double GJPLEPH::GetConstant(GString const_name) {
  return ephHeader.constant[const_name];
}

//
/*****************************************************************************
**                     interp(buf,t,ncf,ncm,na,ifl,pv)                      **
******************************************************************************
**                                                                          **
**    this subroutine differentiates and interpolates a                     **
**    set of chebyshev coefficients to give position and velocity           **
**                                                                          **
**    calling sequence parameters:                                          **
**                                                                          **
**      input:                                                              **
**                                                                          **
**        buf   1st location of array of d.p. chebyshev coefficients        **
**              of position                                                 **
**                                                                          **
**          t   t[0] is double fractional time in interval covered by       **
**              coefficients at which interpolation is wanted               **
**              (0 <= t[0] <= 1).  t[1] is dp length of whole               **
**              interval in input time units.                               **
**                                                                          **
**        ncf   # of coefficients per component                             **
**                                                                          **
**        ncm   # of components per set of coefficients                     **
**                                                                          **
**         na   # of sets of coefficients in full array                     **
**              (i.e., # of sub-intervals in full interval)                 **
**                                                                          **
**         ifl  integer flag: =1 for positions only                         **
**                            =2 for pos and vel                            **
**                                                                          **
**                                                                          **
**      output:                                                             **
**                                                                          **
**        pv   interpolated quantities requested.  dimension                **
**              expected is pv(ncm,ifl), dp.                                **
**                                                                          **
*****************************************************************************/
void GJPLEPH::interp(double coef[], double t[2], int ncf, int ncm, int na,
                     int ifl, double posvel[6]) {
  static double pc[18] = {0.0}, vc[18] = {0.0};  //ÕâÀïÎªÊ²Ã´Ö»ÓÐ18¸öÔªËØ????
  static int np = 2, nv = 3, first = 1;
  static double twot = 0.0;
  double dna, dt1, temp, tc, vfac, temp1;
  int l, i, j;

  if (first) { /* initialize static vectors when called first time */
    pc[0] = 1.0;
    pc[1] = 0.0;
    vc[1] = 1.0;
    first = 0;
  }

  /*  entry point. get correct sub-interval number for this set
   of coefficients and then get normalized chebyshev time
   within that subinterval.                                             */

  dna = (double)na;
  modf(t[0], &dt1);
  temp = dna * t[0];
  l = (int)(temp - dt1);

  /*  tc is the normalized chebyshev time (-1 <= tc <= 1)    */

  tc = 2.0 * (modf(temp, &temp1) + dt1) - 1.0;

  /*  check to see whether chebyshev time has changed,
   and compute new polynomial values if it has.
   (the element pc[1] is the value of t1[tc] and hence
   contains the value of tc on the previous call.)     */

  if (tc != pc[1]) {
    np = 2;
    nv = 3;
    pc[1] = tc;
    twot = tc + tc;
  }

  /*  be sure that at least 'ncf' polynomials have been evaluated
   and are stored in the array 'pc'.    */

  if (np < ncf) {
    for (i = np; i < ncf; ++i) {
      pc[i] = twot * pc[i - 1] - pc[i - 2];
    }

    np = ncf;
  }

  /*  interpolate to get position for each component  */

  for (i = 0; i < ncm; ++i) /* ncm is a number of coordinates */
  {
    posvel[i] = 0.0;
    for (j = ncf - 1; j >= 0; --j) {
      posvel[i] = posvel[i] + pc[j] * coef[j + i * ncf + l * ncf * ncm];
    }
  }

  if (ifl <= 1) return;

  /*  if velocity interpolation is wanted, be sure enough
   derivative polynomials have been generated and stored.    */
  vfac = (dna + dna) / t[1];
  vc[2] = twot + twot;
  if (nv < ncf) {
    for (i = nv; i < ncf; ++i) {
      vc[i] = twot * vc[i - 1] + pc[i - 1] + pc[i - 1] - vc[i - 2];
    }
    nv = ncf;
  }

  /*  interpolate to get velocity for each component    */
  for (i = 0; i < ncm; ++i) {
    posvel[i + ncm] = 0.0;
    for (j = ncf - 1; j > 0; --j) {
      posvel[i + ncm] =
          posvel[i + ncm] + vc[j] * coef[j + i * ncf + l * ncf * ncm];
    }
    posvel[i + ncm] = posvel[i + ncm] * vfac;
  }
  return;
}

/****************************************************************************
****                       split(tt,fr)                                  ****
*****************************************************************************
****  this subroutine breaks a d.p. number into a d.p. integer           ****
****  and a d.p. fractional part.                                        ****
****                                                                     ****
****  calling sequence parameters:                                       ****
****                                                                     ****
****    tt = d.p. input number                                           ****
****                                                                     ****
****    fr = d.p. 2-word output array.                                   ****
****         fr(1) contains integer part                                 ****
****         fr(2) contains fractional part                              ****
****                                                                     ****
****         for negative input numbers, fr(1) contains the next         ****
****         more negative integer; fr(2) contains a positive fraction.  ****
****************************************************************************/
void GJPLEPH::datasplit(double tt, double fr[2]) {
  /*  main entry -- get integer and fractional parts  */
  fr[1] = modf(tt, &fr[0]);
  if (tt >= 0.0 || fr[1] == 0.0) return;
  /*  make adjustments for negative input number   */
  fr[0] = fr[0] - 1.0;
  fr[1] = fr[1] + 1.0;
  return;
}
//

}  // end of namespace gfc
