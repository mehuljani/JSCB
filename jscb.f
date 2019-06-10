	implicit real*8(a-h,o-z)
c	refining the clusters by using positional information/operon
        parameter(ih0=21,ih=84,ifr=3)
	parameter(naa=20)
	parameter(llr=3,mm=4)
	parameter(iigsize=12000000,nngene=12000)
	dimension lend(nngene),irend(nngene)
        dimension lenda(nngene),irenda(nngene)
	dimension iy(iigsize),lena(nngene)
	dimension icsize(nngene)
	dimension iclus(nngene,nngene)
	dimension k(ifr),kaa(nngene)
	dimension rpos(nngene,nngene),irclus(nngene)
	character sign(nngene),signa(nngene),gnome*12000000
    	character*100 genome_file, annotation_file
	character junk
        common/pvt1/ncluscodfreq(nngene,ih),
     .  nclusamfreq(nngene,ih),icluslen(nngene),len(nngene),
     .  ncodonfreq(nngene,ih),naminofreq(nngene,ih)
	common/pvt2/iaa(20,6),iab(20)

c        write(*,*)"Enter genome file name, annotation file name"
c	read(*,*)genome_file,annotation_file

c	nnc=IARGC()
c	call GETARG(1,genome_file)
c        call GETARG(2, annotation_file)

        open(1,file='JSCB_coord')
        open(2,file='JSCBinput')
        open(31,file='JSCB_output.clus')
        open(41,file='JSCB_output.gi')

	nn=0
	do i=1,nngene
	read(1,*,iostat=kk)junk
	if(kk/=0)exit
	nn=nn+1
	enddo
	ngene=nn
	rewind 1

        read(2,*)gnome
	igsize=LEN_TRIM(ADJUSTL(gnome))
	rewind 2

	clusthres=1.d0
	thres2=0.d0
	thres3=0.995d0
        mingenelen=3

        iaa(1,1)=42
        iaa(1,2)=43
        iaa(2,1)=41
        iaa(2,2)=44
        iaa(2,3)=58
        iaa(2,4)=59
        iaa(2,5)=57
        iaa(2,6)=60
        iaa(3,1)=26
        iaa(3,2)=27
        iaa(3,3)=25
        iaa(4,1)=28
        iaa(5,1)=74
        iaa(5,2)=75
        iaa(5,3)=73
        iaa(5,4)=76
        iaa(6,1)=46
        iaa(6,2)=47
        iaa(6,3)=45
        iaa(6,4)=48
        iaa(6,5)=34
        iaa(6,6)=35
        iaa(7,1)=62
        iaa(7,2)=63
        iaa(7,3)=61
        iaa(7,4)=64
        iaa(8,1)=30
        iaa(8,2)=31
        iaa(8,3)=29
        iaa(8,4)=32
        iaa(9,1)=78
        iaa(9,2)=79
        iaa(9,3)=77
        iaa(9,4)=80
        iaa(10,1)=38
        iaa(10,2)=39
        iaa(11,1)=54
        iaa(11,2)=55
        iaa(12,1)=53
        iaa(12,2)=56
        iaa(13,1)=22
        iaa(13,2)=23
        iaa(14,1)=21
        iaa(14,2)=24
        iaa(15,1)=70
        iaa(15,2)=71
        iaa(16,1)=69
        iaa(16,2)=72
        iaa(17,1)=50
        iaa(17,2)=51
        iaa(18,1)=52
        iaa(19,1)=66
        iaa(19,2)=67
        iaa(19,3)=65
        iaa(19,4)=68
        iaa(19,5)=33
        iaa(19,6)=36
        iaa(20,1)=82
        iaa(20,2)=83
        iaa(20,3)=81
        iaa(20,4)=84
                                                                                                                             
        iab(1)=2
        iab(2)=6
        iab(3)=3
        iab(4)=1
        iab(5)=4
        iab(6)=6
        iab(7)=4
        iab(8)=4
        iab(9)=4
        iab(10)=2
        iab(11)=2
        iab(12)=2
        iab(13)=2
        iab(14)=2
        iab(15)=2
        iab(16)=2
        iab(17)=2
        iab(18)=1
        iab(19)=6
        iab(20)=4

        ngenea=0
        do ks=1,ngene
        read(1,*)signa(ks),lenda(ks),irenda(ks)
        lena(ks)=irenda(ks)-lenda(ks)-2
	if(mod(lena(ks),3).ne.0)then
	print*,ks,lena(ks)
	pause
	endif
        if(lena(ks)-3.gt.mingenelen)then
        ngenea=ngenea+1
        sign(ngenea)=signa(ks)
        lend(ngenea)=lenda(ks)
        irend(ngenea)=irenda(ks)
        len(ngenea)=lena(ks)
        endif
        enddo
                                                                                
        thres=.8d0

	do i=1,ngenea
        do j=ih0,ih
        ncodonfreq(i,j)=1
        naminofreq(i,j)=0
        enddo
        enddo
	
	nnp=0
	do ii=1,igsize
        if(gnome(ii:ii).eq.'A')then
        iy(ii)=1
        elseif(gnome(ii:ii).eq.'T')then
        iy(ii)=2
        elseif(gnome(ii:ii).eq.'C')then
        iy(ii)=3
        elseif(gnome(ii:ii).eq.'G')then
        iy(ii)=4
        else
        nnp=nnp+1
c        print*,ii,gnome(ii:ii),iy(ii)
        pause 'problem with genome file'
        endif
        enddo
c        print*,nnp
c        pause

        do j=1,ngenea
        len(j)=irend(j)-lend(j)-2
        if(sign(j).eq.'+')then
        do ii=1,(irend(j)-lend(j)-2)/3
        do ir=1,ifr
        k(ir)=iy(lend(j)+ii*ifr-(ifr-ir)-1)
        enddo
        lx=0
        do ir=1,ifr
        lx=lx+mm**(llr-ir)*k(ir)
        enddo
        ncodonfreq(j,lx)=ncodonfreq(j,lx)+1
        enddo

        do ii=ih0,ih
        do nff=1,naa
        do nff2=1,iab(nff)
        if(iaa(nff,nff2).eq.ii)then
        do nff3=1,iab(nff)
        nff4=iaa(nff,nff3)
        naminofreq(j,ii)=naminofreq(j,ii)+ncodonfreq(j,nff4)
        enddo
        goto 764
        endif
        enddo
        enddo
764     enddo

        elseif(sign(j).eq.'-')then
        do ii=1,(irend(j)-lend(j)-2)/3
        do ir=1,ifr
        if(iy(irend(j)-ii*ifr+(ifr-ir)+1).eq.1)then
        k(ir)=2
        elseif(iy(irend(j)-ii*ifr+(ifr-ir)+1).eq.2)then
        k(ir)=1
        elseif(iy(irend(j)-ii*ifr+(ifr-ir)+1).eq.3)then
        k(ir)=4
        elseif(iy(irend(j)-ii*ifr+(ifr-ir)+1).eq.4)then
        k(ir)=3
        endif
        enddo
        lx=0
        do ir=1,ifr
        lx=lx+mm**(llr-ir)*k(ir)
        enddo
        ncodonfreq(j,lx)=ncodonfreq(j,lx)+1
        enddo
        do ii=ih0,ih
        do nff=1,naa
        do nff2=1,iab(nff)
        if(iaa(nff,nff2).eq.ii)then
        do nff3=1,iab(nff)
        nff4=iaa(nff,nff3)
        naminofreq(j,ii)=naminofreq(j,ii)+ncodonfreq(j,nff4)
        enddo
        goto 765
        endif
        enddo
        enddo
765     enddo
        endif
        enddo

        numgene=1
        nclus=1

52      icsize(nclus)=1
        iclus(nclus,icsize(nclus))=numgene

        do nff2=ih0,ih
        ncluscodfreq(nclus,nff2)=0
        nclusamfreq(nclus,nff2)=0
        enddo

        do nff2=ih0,ih
        do nff=1,icsize(nclus)
        ncluscodfreq(nclus,nff2)=ncluscodfreq(nclus,nff2)
     .  +ncodonfreq(iclus(nclus,nff),nff2)
        enddo
        enddo

        do nff2=ih0,ih
        do nff=1,icsize(nclus)
        nclusamfreq(nclus,nff2)=nclusamfreq(nclus,nff2)
     .  +naminofreq(iclus(nclus,nff),nff2)
        enddo
        enddo
        if(numgene.eq.ngenea)then
        nclus0=nclus
        write(99,*)'nclus0',nclus0
        goto 50
        endif

51      icluslen(nclus)=0
        do ii=1,naa
        icluslen(nclus)=icluslen(nclus)+nclusamfreq(nclus,iaa(ii,1))
        enddo
        numgene=numgene+1
        lengtot=icluslen(nclus)+len(numgene)/3
        rnn=lengtot
        ddiv=divgn2(numgene,nclus,lengtot)
        pbc=sig(rnn,abs(ddiv))
        if(pbc.lt.thres)then

        do nff2=ih0,ih
        ncluscodfreq(nclus,nff2)=ncluscodfreq(nclus,nff2)+
     .  ncodonfreq(numgene,nff2)
        nclusamfreq(nclus,nff2)=nclusamfreq(nclus,nff2)+
     .  naminofreq(numgene,nff2)
        enddo

        icsize(nclus)=icsize(nclus)+1
        iclus(nclus,icsize(nclus))=numgene
        if(numgene.eq.ngenea)then
        nclus0=nclus
        write(99,*)'nclus0',nclus0
        goto 50
        endif
        goto 51

        else

        nclus=nclus+1
        goto 52
        endif

50	do i=1,nclus
        icluslen(i)=0
	enddo
	do i=1,nclus
	do ii=1,naa
	icluslen(i)=icluslen(i)+nclusamfreq(i,iaa(ii,1))
	enddo
	enddo

        do kkn1=1,nclus-1
        do kkn2=kkn1+1,nclus
	lengtot=icluslen(kkn1)+icluslen(kkn2)
        djs=divgn(kkn1,kkn2,lengtot)
	rnn=icluslen(kkn1)+icluslen(kkn2)
	pbc=sig(rnn,djs)
	if(pbc.lt.thres3)then
	do nff2=ih0,ih
	ncluscodfreq(kkn1,nff2)=ncluscodfreq(kkn1,nff2)+
     .  ncluscodfreq(kkn2,nff2)
	nclusamfreq(kkn1,nff2)=nclusamfreq(kkn1,nff2)+
     .  nclusamfreq(kkn2,nff2)
	enddo

	do i=1,icsize(kkn2)
	iclus(kkn1,icsize(kkn1)+i)=iclus(kkn2,i)
	enddo	
	
	icsize(kkn1)=icsize(kkn1)+icsize(kkn2)

	do i=kkn2,nclus-1
	do j=1,icsize(i+1)
	iclus(i,j)=iclus(i+1,j)
	enddo
	enddo
	
	do i=kkn2,nclus-1
	icsize(i)=icsize(i+1)
	do j=ih0,ih
	ncluscodfreq(i,j)=ncluscodfreq(i+1,j)
	nclusamfreq(i,j)=nclusamfreq(i+1,j)
	enddo
	enddo
	nclus=nclus-1
	goto 50
	endif
	enddo
        enddo

59      do i=1,nclus
        do j=1,nclus
        rpos(i,j)=0.d0
        enddo
        enddo

        do i=1,nclus
        do j=1,icsize(i)
        do ii=1,nclus
        do jj=1,icsize(ii)
        if(iclus(i,j)-1.eq.iclus(ii,jj))then
        do iii=1,nclus
        do jjj=1,icsize(iii)
        if(iclus(i,j)+1.eq.iclus(iii,jjj))then
        if(ii.eq.iii)then
        rpos(i,ii)=rpos(i,ii)+1.d0
        goto 779
        endif
        endif
        enddo
        enddo
        endif
        enddo
        enddo
779     enddo
        enddo
                                                                                

996     do i=1,nclus-1
        do j=i+1,nclus
c        rp=(rpos(i,j)/(icsize(i)-rpos(i,i))+
c     .  rpos(j,i)/(icsize(j)-rpos(j,j)))/2.d0
        rp=(rpos(i,j)/icsize(i)+
     .  rpos(j,i)/icsize(j))/2.d0
	if(rp.gt.clusthres)then

        do nff2=ih0,ih
        ncluscodfreq(i,nff2)=ncluscodfreq(i,nff2)+
     .  ncluscodfreq(j,nff2)
        nclusamfreq(i,nff2)=nclusamfreq(i,nff2)+
     .  nclusamfreq(j,nff2)
        enddo
                                                                                
	do m=1,icsize(j)
	iclus(i,icsize(i)+m)=iclus(j,m)
	enddo
	icsize(i)=icsize(i)+icsize(j)

        do iu=j,nclus-1
        do ju=1,icsize(iu+1)
        iclus(iu,ju)=iclus(iu+1,ju)
        enddo
        enddo

        do iu=j,nclus-1
        icsize(iu)=icsize(iu+1)
	do ju=ih0,ih
        ncluscodfreq(iu,ju)=ncluscodfreq(iu+1,ju)
        nclusamfreq(iu,ju)=nclusamfreq(iu+1,ju)
        enddo
	enddo

	nclus=nclus-1

        do mi=1,nclus
        do mj=1,nclus
        rpos(mi,mj)=0.d0
        enddo
        enddo
                                                                                
        do mi=1,nclus
        do mj=1,icsize(mi)
        do ii=1,nclus
        do jj=1,icsize(ii)
        if(iclus(mi,mj)-1.eq.iclus(ii,jj))then
        do iii=1,nclus
        do jjj=1,icsize(iii)
        if(iclus(mi,mj)+1.eq.iclus(iii,jjj))then
        if(ii.eq.iii)then
        rpos(mi,ii)=rpos(mi,ii)+1.d0
        goto 778
        endif
        endif
        enddo
        enddo
        endif
        enddo
        enddo
778     enddo
        enddo

	goto 996
	endif
        enddo
        enddo

	do i=1,nclus
        icluslen(i)=0
        enddo
        do i=1,nclus
        do ii=1,naa
        icluslen(i)=icluslen(i)+nclusamfreq(i,iaa(ii,1))
        enddo
        enddo

	nnh=0
	nnm=0

        do i=1,nclus
	kw=0
        do j=1,icsize(i)
	jfk=0
        do ii=1,nclus
        do jj=1,icsize(ii)
	if(iclus(i,j).eq.1.or.iclus(i,j).eq.ngenea)goto 299
        if(iclus(i,j)-1.eq.iclus(ii,jj))then
	mgg1=ii
	jfk=jfk+1
        elseif(iclus(i,j)+1.eq.iclus(ii,jj))then
        mgg2=ii
	jfk=jfk+1
	endif
	if(jfk.eq.2)goto 29
        enddo
        enddo
29      if(mgg1.eq.mgg2.and.mgg1.ne.i)then
	nnh=nnh+1
        lengtot=len(iclus(i,j))/3+icluslen(mgg1)
        gdjs=divgn2(iclus(i,j),mgg1,lengtot)
	rnn=lengtot
        pbc=sig(rnn,gdjs)
c        if(pbc.lt.thres)then
        if(pbc.lt.thres2)then
	nnm=nnm+1
	kw=kw+1
	icsize(mgg1)=icsize(mgg1)+1
	iclus(mgg1,icsize(mgg1))=iclus(i,j)
	icluslen(mgg1)=icluslen(mgg1)+len(iclus(i,j))/3
        do nff2=ih0,ih
        ncluscodfreq(mgg1,nff2)=ncluscodfreq(mgg1,nff2)+
     .  ncodonfreq(iclus(i,j),nff2)
        nclusamfreq(mgg1,nff2)=nclusamfreq(mgg1,nff2)+
     .  naminofreq(iclus(i,j),nff2)
        enddo

	kaa(kw)=j
	endif
	endif
299	enddo
	kdd=0
	do ii=1,icsize(i)
	do jj=1,kw
	if(ii.eq.kaa(jj))goto 399
	enddo
	kdd=kdd+1
	irclus(kdd)=iclus(i,ii)
399	enddo
	icsize(i)=kdd
	do ii=1,kdd
	iclus(i,ii)=irclus(ii)
	enddo

        do nff2=ih0,ih
        ncluscodfreq(i,nff2)=0
        nclusamfreq(i,nff2)=0
        enddo
                                                                                
        do nff2=ih0,ih
        do nff=1,icsize(i)
        ncluscodfreq(i,nff2)=ncluscodfreq(i,nff2)
     .  +ncodonfreq(iclus(i,nff),nff2)
        enddo
        enddo
                                                                                
        do nff2=ih0,ih
        do nff=1,icsize(i)
        nclusamfreq(i,nff2)=nclusamfreq(i,nff2)
     .  +naminofreq(iclus(i,nff),nff2)
        enddo
        enddo

        icluslen(i)=0
        do ii=1,naa
        icluslen(i)=icluslen(i)+nclusamfreq(i,iaa(ii,1))
        enddo

	enddo

        do i=1,nclus
        do j=1,icsize(i)
        do m=1,ngene
        if(lenda(m).eq.lend(iclus(i,j)).and.irenda(m)
     .  .eq.irend(iclus(i,j)))then
        write(31,*)i,icsize(i),m,lena(m)
        endif
        enddo
        enddo
        enddo

        do m=1,ngene
        do i=1,nclus
        do j=1,icsize(i)
        if(lenda(m).eq.lend(iclus(i,j)).and.irenda(m)
     .  .eq.irend(iclus(i,j)))then
        write(41,*)m,i,icsize(i),lena(m)
        endif
        enddo
        enddo
        enddo

	end



	FUNCTION rentrotot(id1,id2,lengtot)
	implicit real*8(a-h,o-z)
	parameter(nngene=12000)
	parameter(ih0=21,ih=84,naa=20)
	common/pvt1/ncluscodfreq(nngene,ih),
     .  nclusamfreq(nngene,ih),icluslen(nngene),len(nngene),
     .  ncodonfreq(nngene,ih),naminofreq(nngene,ih)
	common/pvt2/iaa(20,6),iab(20)
	rentrotot=0.d0
        do nff=1,naa
        rentroi=0.d0
        do nff2=1,iab(nff)
        pp=(ncluscodfreq(id1,iaa(nff,nff2))+
     .  ncluscodfreq(id2,iaa(nff,nff2)))/dfloat
     .  (nclusamfreq(id1,iaa(nff,1))+nclusamfreq(id2,iaa(nff,1)))
        rentroi=rentroi+pp*log(pp)
        enddo
        pp0=(nclusamfreq(id1,iaa(nff,1))+nclusamfreq(id2,iaa(nff,1)))/
     .  dfloat(lengtot)
        rentroi=rentroi*pp0
        rentrotot=rentrotot+rentroi
        enddo
        rentrotot=-rentrotot
        return
        end

        FUNCTION rentrotot2(id1,id2,lengtot)
        implicit real*8(a-h,o-z)
        parameter(nngene=12000)
        parameter(ih0=21,ih=84,naa=20)
        common/pvt1/ncluscodfreq(nngene,ih),
     .  nclusamfreq(nngene,ih),icluslen(nngene),len(nngene),
     .  ncodonfreq(nngene,ih),naminofreq(nngene,ih)
        common/pvt2/iaa(20,6),iab(20)
        rentrotot2=0.d0
        do nff=1,naa
        rentroi=0.d0
        do nff2=1,iab(nff)
        pp=(ncodonfreq(id1,iaa(nff,nff2))+
     .  ncluscodfreq(id2,iaa(nff,nff2)))/dfloat
     .  (naminofreq(id1,iaa(nff,1))+nclusamfreq(id2,iaa(nff,1)))
        rentroi=rentroi+pp*log(pp)
        enddo
        pp0=(naminofreq(id1,iaa(nff,1))+nclusamfreq(id2,iaa(nff,1)))/
     .  dfloat(lengtot)
        rentroi=rentroi*pp0
        rentrotot2=rentrotot2+rentroi
        enddo
        rentrotot2=-rentrotot2
        return
        end



	FUNCTION rentro(id)
	implicit real*8(a-h,o-z)
	parameter(nngene=12000)
	parameter(ih0=21,ih=84,naa=20)
	common/pvt1/ncluscodfreq(nngene,ih),
     .  nclusamfreq(nngene,ih),icluslen(nngene),len(nngene),
     .  ncodonfreq(nngene,ih),naminofreq(nngene,ih)
	common/pvt2/iaa(20,6),iab(20)
        rentro=0.d0
	do nff=1,naa
	rentroi=0.d0
        do nff2=1,iab(nff)
	pp=ncluscodfreq(id,iaa(nff,nff2))/
     .  dfloat(nclusamfreq(id,iaa(nff,nff2)))
	rentroi=rentroi+pp*log(pp)
	enddo
       rentroi=rentroi*(nclusamfreq(id,iaa(nff,1))/dfloat(icluslen(id)))
	rentro=rentro+rentroi
	enddo
	rentro=-rentro
        return
        end


        FUNCTION rentro2(id)
        implicit real*8(a-h,o-z)
        parameter(nngene=12000)
        parameter(ih0=21,ih=84,naa=20)
        common/pvt1/ncluscodfreq(nngene,ih),
     .  nclusamfreq(nngene,ih),icluslen(nngene),len(nngene),
     .  ncodonfreq(nngene,ih),naminofreq(nngene,ih)
        common/pvt2/iaa(20,6),iab(20)
        rentro2=0.d0
        do nff=1,naa
        rentroi=0.d0
        do nff2=1,iab(nff)
        pp=ncodonfreq(id,iaa(nff,nff2))/
     .  dfloat(naminofreq(id,iaa(nff,nff2)))
        rentroi=rentroi+pp*log(pp)
        enddo
       rentroi=rentroi*(naminofreq(id,iaa(nff,1))/dfloat(len(id)/3))
        rentro2=rentro2+rentroi
        enddo
        rentro2=-rentro2
        return
        end


        FUNCTION divgn(i,j,lengtot)
	implicit real*8(a-h,o-z)
	parameter(nngene=12000)
	parameter(ih0=21,ih=84)
	common/pvt1/ncluscodfreq(nngene,ih),
     .  nclusamfreq(nngene,ih),icluslen(nngene),len(nngene),
     .  ncodonfreq(nngene,ih),naminofreq(nngene,ih)
        hhh=2.0d0
        xxx=1/log(hhh)
        divgn=xxx*(rentrotot(i,j,lengtot)-
     .  (icluslen(i)*rentro(i)+icluslen(j)*rentro(j))/
     .  dfloat(lengtot))
        return
        end


        FUNCTION divgn2(i,j,lengtot)
        implicit real*8(a-h,o-z)
        parameter(nngene=12000)
        parameter(ih0=21,ih=84)
        common/pvt1/ncluscodfreq(nngene,ih),
     .  nclusamfreq(nngene,ih),icluslen(nngene),len(nngene),
     .  ncodonfreq(nngene,ih),naminofreq(nngene,ih)
        hhh=2.0d0
        xxx=1/log(hhh)
        divgn2=xxx*(rentrotot2(i,j,lengtot)-
     .  ((len(i)/3)*rentro2(i)+icluslen(j)*rentro(j))/
     .  dfloat(lengtot))
        return
        end


        FUNCTION sig(rnn,djsmx)
	implicit real*8(a-h,o-z)
	snn=rnn*log(2.d0)*djsmx
	sig=gammp(30.d0,snn)
        return
        end

	

	 FUNCTION gammp(a,x)
      REAL*8 a,gammp,x
C     USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.) then
c	write(101,*)a,x	
	pause 'bad arguments in gammp'
	endif
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
       call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      END


      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=9000,EPS=3.e-7)
C    USES gammln
      INTEGER n
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END


      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=9000,EPS=3.e-7,FPMIN=1.e-30)
C     USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
	 if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END



      FUNCTION gammln(xx)
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282623310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END


      FUNCTION ran3a(nth,id,idum)
      INTEGER idum,id,nth
      INTEGER MBIG,MSEED,MZ,ik,nthres
C     REAL MBIG,MSEED,MZ
      REAL ran3a,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      PARAMETER (nthres=20,ik=20)
      INTEGER i,iff(nthres,ik),ii,inext(nthres,ik),inextp(nthres,ik),k
      INTEGER mj,mk,ma(nthres,ik,55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff / 400*0 /
      if(idum.lt.0.or.iff(nth,id).eq.0)then
        iff(nth,id)=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(nth,id,55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(nth,id,ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
        mj=ma(nth,id,ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(nth,id,i)=ma(nth,id,i)-ma(nth,id,1+mod(i+30,55))
            if(ma(nth,id,i).lt.MZ)ma(nth,id,i)=ma(nth,id,i)+MBIG
12        continue
13      continue
        inext(nth,id)=0
        inextp(nth,id)=31
        idum=1
      endif
      inext(nth,id)=inext(nth,id)+1
      if(inext(nth,id).eq.56)inext(nth,id)=1
      inextp(nth,id)=inextp(nth,id)+1
      if(inextp(nth,id).eq.56)inextp(nth,id)=1
      mj=ma(nth,id,inext(nth,id))-ma(nth,id,inextp(nth,id))
      if(mj.lt.MZ)mj=mj+MBIG
      ma(nth,id,inext(nth,id))=mj
      ran3a=mj*FAC
      return
      END

      FUNCTION ran3b(nth,id,idum)
      INTEGER nth,idum,id
      INTEGER MBIG,MSEED,MZ,ik,nthres
C     REAL MBIG,MSEED,MZ
      REAL ran3b,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      PARAMETER (nthres=20,ik=20)
      INTEGER i,iff(nthres,ik),ii,inext(nthres,ik),inextp(nthres,ik),k
      INTEGER mj,mk,ma(nthres,ik,55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff / 400*0 /
      if(idum.lt.0.or.iff(nth,id).eq.0)then
        iff(nth,id)=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(nth,id,55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(nth,id,ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
        mj=ma(nth,id,ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(nth,id,i)=ma(nth,id,i)-ma(nth,id,1+mod(i+30,55))
            if(ma(nth,id,i).lt.MZ)ma(nth,id,i)=ma(nth,id,i)+MBIG
12        continue
13      continue
        inext(nth,id)=0
        inextp(nth,id)=31
        idum=1
      endif
      inext(nth,id)=inext(nth,id)+1
      if(inext(nth,id).eq.56)inext(nth,id)=1
      inextp(nth,id)=inextp(nth,id)+1
      if(inextp(nth,id).eq.56)inextp(nth,id)=1
      mj=ma(nth,id,inext(nth,id))-ma(nth,id,inextp(nth,id))
      if(mj.lt.MZ)mj=mj+MBIG
      ma(nth,id,inext(nth,id))=mj
      ran3b=mj*FAC
      return
      END


      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=500)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

