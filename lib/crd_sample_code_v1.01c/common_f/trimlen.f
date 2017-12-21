      integer function trimlen(array)
C  find the length of an array to the first blank (or end of array)
C
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      character array*(*)

      length= len(array)
      trimlen= length
      do 10 i=1,length
        lenbk= length+1-i
        if (array(lenbk:lenbk).eq.' ') go to 10
        if (array(lenbk:lenbk).eq.char(0)) go to 10
CC        if (array(lenbk:lenbk).eq.' ' .or.
CC     .          array(lenbk:lenbk).eq.char(0)) go to 10
          trimlen= lenbk
          go to 90
10    continue
      trimlen= 0
90    return
      end

