      program fort_ratio
      real yes(15,91,11),noo(15,91,11),ratio(15,91,11),jxdiff
      real altyes(21,91,16),altnoo(21,91,16),altratio(21,91,16)

      open(unit=101, file="long_slice_000yesboth.txt")
      open(unit=102, file="long_slice_090yesboth.txt")
      open(unit=103, file="long_slice_180yesboth.txt")
      open(unit=104, file="long_slice_270yesboth.txt")
      open(unit=105, file="alt_slice_k-10yes.txt")

      open(unit=106, file="long_slice_000nooboth.txt")
      open(unit=107, file="long_slice_090nooboth.txt")
      open(unit=108, file="long_slice_180nooboth.txt")
      open(unit=109, file="long_slice_270nooboth.txt")
      open(unit=110, file="alt_slice_k-10noo.txt")

      open(unit=111, file="long_slice_000rat.txt")
      open(unit=112, file="long_slice_090rat.txt")
      open(unit=113, file="long_slice_180rat.txt")
      open(unit=114, file="long_slice_270rat.txt")
      open(unit=115, file="alt_slice_k-10rat.txt")

      do I=1,58
      do J=1,15
         read(101,*) yes(J,I,1),yes(J,I,2),yes(J,I,3),yes(J,I,4),
     1yes(J,I,5),yes(J,I,6),yes(J,I,7),yes(J,I,8),yes(J,I,9),
     2yes(J,I,10),yes(J,I,11)    
      enddo
c      read(101,*)dummy
      enddo

      do I=1,58
      do J=1,15
         read(106,*) noo(J,I,1),noo(J,I,2),noo(J,I,3),noo(J,I,4),
     1noo(J,I,5),noo(J,I,6),noo(J,I,7),noo(J,I,8),noo(J,I,9),
     2noo(J,I,10),noo(J,I,11)    
      enddo
      enddo

      do I=1,58
      do J=1,15
      write(111,*) noo(J,I,1),noo(J,I,2),noo(J,I,3),
     1noo(J,I,4),yes(J,I,5)/noo(J,I,5),yes(J,I,6)/noo(J,I,6),
     2yes(J,I,7)/noo(J,I,7),yes(J,I,8)/noo(J,I,8),
     3abs(yes(J,I,9)-noo(J,I,9)),abs(yes(J,I,10)-noo(J,I,10)),
     4yes(J,I,11)/noo(J,I,11) 
      enddo
      write(111,*)" "
      enddo

      do I=1,91
      do J=1,21
         read(105,*) altyes(J,I,1),altyes(J,I,2),altyes(J,I,3),
     1altyes(J,I,4),altyes(J,I,5),altyes(J,I,6),altyes(J,I,7),
     2altyes(J,I,8),altyes(J,I,9),altyes(J,I,10),altyes(J,I,11),
     3altyes(J,I,12),altyes(J,I,13),altyes(J,I,14),altyes(J,I,15),
     4altyes(J,I,16)    
      enddo
      enddo

      do I=1,91
      do J=1,21
         read(110,*) altnoo(J,I,1),altnoo(J,I,2),altnoo(J,I,3),
     1altnoo(J,I,4),altnoo(J,I,5),altnoo(J,I,6),altnoo(J,I,7),
     2altnoo(J,I,8),altnoo(J,I,9),altnoo(J,I,10),altnoo(J,I,11),
     3altnoo(J,I,12),altnoo(J,I,13),altnoo(J,I,14),altnoo(J,I,15),
     4altnoo(J,I,16)    
      enddo
      enddo

      do I=1,91
      do J=1,21
         if ((altyes(J,I,12) .gt. 0) .and. (altnoo(J,I,12) .gt. 0))then
            jxdif = (altyes(J,I,12)-altnoo(J,I,12)) * 1
      else if ((altyes(J,I,12) .lt. 0).and.(altnoo(J,I,12) .lt. 0))then
            jxdif = (altyes(J,I,12)-altnoo(J,I,12)) * (-1)
         else
            jxdif = altyes(J,I,12)-altnoo(J,I,12)
            
         endif
               
         write(115,*) altyes(J,I,1),altyes(J,I,2),altyes(J,I,3),
     1altyes(J,I,4),altyes(J,I,5),altyes(J,I,6),
     2altyes(J,I,7)/altnoo(J,I,7),altyes(J,I,8)/altnoo(J,I,8),
     3altyes(J,I,9)/altnoo(J,I,9),altyes(J,I,10)/altnoo(J,I,10),
     4altyes(J,I,11)/altnoo(J,I,11),(altyes(J,I,12)-altnoo(J,I,12)),
     5(altyes(J,I,13)-altnoo(J,I,13)),altyes(J,I,14)/altnoo(J,I,14),
     6altyes(J,I,15)/altnoo(J,I,15),altyes(J,I,16)/altnoo(J,I,16)
      enddo
      write(115,*)" "
      enddo




      end program
