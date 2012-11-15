module dannie
character(20) npr !название проекта
integer(4) NV !количество веществ
character(1000) NAZ(10) !названия файлов !можно увеличить до 100
character(1000) NAZF(10) !названия файлов
real(8) RoNV,RoNVML !Молекулярная плотность
real(8) V !Элементарный объем
integer(4) N !Количество молекул в объеме
integer(4) storona !количество узлов кристаллической решетки
integer(1) tipResh !Тип решетки 1-объемноцентриованная 2-гранецентрированная
real(8) DlSt ! длинна элементарного кубика
real(8) KonSt !координата конца стороны
integer(4) NmMax !максимальное количество атомов в молекуле
integer(4) totalNa !общее количество атомов
real(8) sproc(100) !состав в процентах
real(8), allocatable:: x(:) !x координата
real(8), allocatable:: y(:) !y координата
real(8), allocatable:: z(:) !z координата
integer(4), allocatable:: kolNv(:) !количество атомов данного типа
!координаты атомов в молекулеl
real(8) xm(100,100) !координата х атома в молекуле
real(8) ym(100,100) !координата у атома в молекуле
real(8) zm(100,100) !координата z атома в молекуле
real(8) rm(100,100) !размеры атомов
character(5) labelm(100,100) !название атома
character(5), allocatable:: labela(:,:)
!виртуальные координаты атома для проверки растояния
real(8) xprov1(100) !координата х атома в молекуле
real(8) yprov1(100) !координата у атома в молекуле
real(8) zprov1(100) !координата z атома в молекуле
!
real(8) bs,be !безразмерный эпсилон и сигма
real(8) xprov2(100) !координата х атома в молекуле
real(8) yprov2(100) !координата у атома в молекуле
real(8) zprov2(100) !координата z атома в молекуле
integer(4) Nm(100) !количество атомов в молекуле
real(8) MM(100)
integer(4),allocatable:: Na(:) !количество центров в конкретной молекуле
integer(4),allocatable:: Ntip(:) !тип (номер) атома
!Координаты относительно центра молекулы (относительные)
real(8),allocatable:: xa(:,:) !координаты x атома
real(8),allocatable:: ya(:,:) !координаты y атома
real(8),allocatable:: za(:,:) !координаты z атома
!Координаты относительно центра молекулы (абсолютные)
real(8),allocatable:: xaa(:,:) !координаты x атома
real(8),allocatable:: yaa(:,:) !координаты y атома
real(8),allocatable:: zaa(:,:) !координаты z атома
integer(4),allocatable:: koltors(:) !количество торсионных связей
real(8),allocatable:: CosTorsU(:,:) !координаты z атома
real(8),allocatable:: TorsU(:,:)    !торсионный угол

real(8) CosTorsU_do(30)
real(8) TorsU_do(30)

real(8):: start !начало отсчета
real(8):: finish !конец отсчета
real(8):: time1
real(8):: time2

real(8) xdo !координата молекулы до перемещения
real(8) ydo !координата молекулы до перемещения
real(8) zdo !координата молекулы до перемещения

real(8),allocatable:: xado(:) !координаты x атома
real(8),allocatable:: yado(:) !координаты y атома
real(8),allocatable:: zado(:) !координаты z атома
!real(8) xado(30) !координата атома до перемещения
!real(8) yado(30) !координата атома до перемещения
!real(8) zado(30) !координата атома до перемещения

real(8),allocatable:: xaado(:) !координаты x атома
real(8),allocatable:: yaado(:) !координаты y атома
real(8),allocatable:: zaado(:) !координаты z атома
!real(8) xaado(30) !координата атома до перемещения
!real(8) yaado(30) !координата атома до перемещения
!real(8) zaado(30) !координата атома до перемещения
real(8) siga(30,30,30,30) !параметры центр-центрового взаимодействия
real(8) epsa(30,30,30,30)
real(8) alf(30,30,30,30)
real(8) kk1(30,30,30,30)
real(8) kk2(30,30,30,30)
!real(8) kk3(30,30,30,30)
!real(8) kk4(30,30,30,30)
real(8) dkk1(30,30,30,30)
real(8) dkk2(30,30,30,30)
real(8) dkk3(30,30,30,30)
real(8) dkk4(30,30,30,30)
real(8) siga6
real(8) bk1(30,30,30,30)
real(8) bk2(30,30,30,30)
real(8) dbk1(30,30,30,30)
real(8) dbk2(30,30,30,30)

character(20) MMname(30)
real(8) rz2,rz6,rz,rz7,expp,expp6
real(8) DlShag !длинна масимального шага
real(8) randn !случайное число для нахождения передвигаемой молекулы
real(8) rastmm !Растояние между двумя молекулами
real(8) epsi(100,100) !эпсилон
real(8) sigma(100,100) !сигма
real(8) alfa(100,100) !альфа букингема
real(8) dopeps(100,100,100,100)
integer(4) m1,m2,a1,a2 !2молекулы и атома
real(8) potenc_do, potenc_posle !потенциал до и после перемещения
real(8) vir_do, vir_posle
integer(4) assept,notassept !количество принятых и не принтых перемещений
real(8) Temp,KoefP,TempK !температура безразмерная, коэффициент поворота
real(8) RadCut !радиус обрезания
real(8) LJatom,LJmol, LJpotenc !межатомный, межмолекулярный, для одной молекулы с окружением
real(8) deltaP,deltaVir !разница потенциала
real(8) sumP, TotEn, en ,ven,vven
integer(4) as, aDl, naDl
real(8) rastaa
real(8) dmx, dmy, dmz !проекции растояния между молекулами
real(8) dax, day,daz !проекции растояния между атомами
real(8) Vir, VirLJatom, VirLJmol !вириал межмолекулярный, межатомный
real(8) TotalVir,vi,sumVir
real(8) dx,dy,dz,d
integer(4) Nout, Nmov, Nchan, Neqv,Ntot, Nnew
integer(4) mov !номер перемещения
real(8) LJtail,LJPtail ,RadCut2 !хвост потенциала
real(8) ChPI, ChPI2 !ПИ
real(8) deltar !для построения гистограммы
real(8),allocatable:: Hist(:) !гистограмма
real(8),allocatable:: Gr(:) !гистограмма
real(8),allocatable:: SumGr(:) !сумма для вириала гистограмма
real(8),allocatable:: VirGr(:) !вириал гистограммы
real(8),allocatable:: Nideal(:) !Идеальное распределение
real(8),allocatable:: Grx(:) !гистограмма
!для атом -атомного распределения
real(8),allocatable:: HistAtom(:,:) !гистограмма
real(8),allocatable:: GrA(:,:) !гистограмма
real(8),allocatable:: SumGrA(:,:) !сумма для вириала гистограмма
real(8),allocatable:: VirGrA(:,:) !вириал гистограммы
real(8),allocatable:: SumGrAN(:,:) !сумма для вириала гистограмма
real(8),allocatable:: VirGrAN(:,:) !вириал гистограммы
real(8),allocatable:: RDF_Block(:,:,:)
!
real(8) GRconst !константа для гистограммы
real(8) rVerh, rNiz, rGrcut
integer(4) HistRaz !размер гистограммы
real(8) davl, enid, sen, svi, sdavl
integer(4) TipCut, TipLat !тип обрезания и считывания начальной структуры
integer(4) anom !номер атома
real(8) sumPout, sumVirout, enout,viout,davlout,senout,sdavlout
integer(4) Neqvil, Mout, Meqvil
real(8) drtail,ptail1,ptail2,etail1,etail2
integer(4) nstep
integer(4) atomk !количество атомов в системе
real(8),allocatable:: oenp(:) !хранятся промежуточные значения энергии
real(8),allocatable:: odavl(:) !хранятся промежуточные значения давления
real(8),allocatable:: ocv(:)
real(8) outsenout, outdavlout,o2enp,o2davl
character(200) workdir
integer(2) ena
character(1000) Nazv
character(1100) Nazv2,Nazv3
real(8) razme, razmd
integer(4) Sample,SampleEqv
integer(4) ptip !выбор типа потенциала
integer(4) Ntors(30)
integer(4) fsta(30,30)
integer(4) seca(30,30)
integer(4) nrota(30,30)
integer(4) rota(30,30,30)
integer(4) srota(30,30,30)
real(8) torskoef(30,30,30)
real(8) torscos(30,30)
real(8) torsEnergy !торсионная энергия
integer(4) NtMax !максимальное количество торсионных поворотов в момлекуле
integer(4) Mcm(30)
real(8) xcm(30,100),ycm(30,100),zcm(30,100)
character(5) labm(30,100) !название атома
integer(4) kolcent !общее количество цетров в атомах всех типов
character(5),allocatable:: cnames(:)
integer(4),allocatable:: nnames(:)
integer(4) nomer, koln, ckol
character(5),allocatable:: labckol(:)
integer(4) kolfrr, inom, jnom
character(100) rdffile,rdffile2
integer(4) Mrdf
integer(4),allocatable:: Ncrdf(:,:)
real(8),allocatable:: Ncsproc(:)
real(8) Ncsum1, Ncsum2
integer(4) sumKNv1, sumKNv2
real(8),allocatable:: y_rast(:,:)
real(8),allocatable:: HistY(:,:)
integer(4),allocatable:: nom_y(:,:)
real(8) potencm,potenc_minus
character(100) yfile
real(8),allocatable:: sum_mu(:)
integer(4),allocatable:: nom_mu(:)
real(8),allocatable:: mu_out(:)
real(8),allocatable:: sum_minus(:,:)
integer(4) Nk_rand,kon_rast
real(8),allocatable:: y_potenc2(:,:)
real(8),allocatable:: sum_potenc2(:,:)
real(8) delta_rast
integer(4) y2_nom
real(8) beze,bezro,bezp,bezt
character(20) Razm_davl,razm_ro,razm_tem,Razm_eN,razm_rast,razm_mu,razm_cv
real(8) bezcv
real(8),allocatable:: total_c_f(:,:)
real(8),allocatable:: AFFT(:)
real(8),allocatable:: BFFT(:)
real(8),allocatable:: WorkA(:)
real(8),allocatable:: WorkB(:)
integer(4) stepen
real(8) deltak
real(8),allocatable:: direct_c_f(:,:)
real(8),allocatable:: cavity_c_f(:,:)
integer(4),allocatable:: mu2_nom(:)
real(8),allocatable:: mu2_sum(:)
real(8),allocatable:: mu2_out(:)
real(8),allocatable:: cavity_out(:,:)
real(8) maxsig
integer(4) kol_m
real(8),allocatable:: x_m(:)
real(8),allocatable:: y_m(:)
real(8),allocatable:: z_m(:)
real(8),allocatable:: xaa_m(:,:)
real(8),allocatable:: yaa_m(:,:)
real(8),allocatable:: zaa_m(:,:)
real(8),allocatable:: HistAtom_m(:,:)
real(8),allocatable:: GrA_m(:,:)
real(8),allocatable:: VirGrA_m(:,:)
real(8),allocatable:: Grx_m(:)
real(8),allocatable:: Nideal_m(:)
real(8),allocatable:: total_c_f_m(:,:)
real(8),allocatable:: indirect_c_f(:,:)
real(8),allocatable:: SumGrA_m(:,:)
real(8),allocatable:: Bridge_f(:,:)
real(8),allocatable:: smooth(:)
real(8) rGrcut_m
integer(4) HistRaz_m, N_m, Mrdf_m, MRDFN
integer(4) calc_rdf,calc_hend
real(8) sec_t
integer(4) s_i,s_k
integer(4) hour_t,min_t, sm_i
character(10) versia
integer(4) bind_kol(30)
integer(4) bind_cent(30,30)
integer(4) bind1(30,30)
integer(4) bind2(30,30)
integer(4) bind1k(30,30,30)
integer(4) bind2k(30,30,30)
real(8) bind_k(30,30)
real(8) bind_fi(30,30)
real(8),allocatable:: bind_u(:,:)     !угол в радианах
real(8) xv1,xv2,yv1,yv2,zv1,zv2
real(8) r1_bind, r2_bind,cos_bind
real(8) bind_en
real(8) bind_u_do(30)
integer(4) prov_var,kol_deg
character(200) bindfile
integer(4),allocatable:: hist_bind(:,:,:)
integer(4),allocatable:: hist_norm(:,:)
real(8),allocatable:: hist_u(:,:,:)
real(8) deg,deg_sh
integer(4) sum_hist
!
integer(4) l_kol(30)
real(8) l_k(30)
real(8) l_r(30)
integer(4) l_n1(30)
integer(4) l_n2(30)
integer(4) l_nk1(30,30)
integer(4) l_nk2(30,30)
!
integer(4) bind_val(30)
integer(4) bind_vn(30,30,3)
real(8) x1valen,x2valent,y1valent,y2valent,z1valent,z2valent
real(8) valcos
integer(4) bond_kol(30)
integer(4) bond_n(30,30,2)
!заряды
integer(4) charg_kol(30)
real(8) ch_x(30,30)            !координаты зарядов в молекуле
real(8) ch_y(30,30)
real(8) ch_z(30,30)
real(8) ch_q(30,30)
real(8),allocatable:: ch_ax(:,:)         !координаты конкретных зарядов
real(8),allocatable:: ch_ay(:,:)
real(8),allocatable:: ch_az(:,:)
real(8),allocatable:: ch_aq(:,:)
real(8),allocatable:: ch_aax(:,:)        !аболютные координаты зарядов
real(8),allocatable:: ch_aay(:,:)
real(8),allocatable:: ch_aaz(:,:)
integer(4) NmCh
integer(4) ch_bind1(30,30)           !кручение зарядов
integer(4) ch_bind2(30,30)
integer(4) ch_bind1k(30,30,30)       !
integer(4) ch_bind2k(30,30,30)
real(8) rb_kc(30,30)                 !граница определения КЧ
real(8) kc(30,30)
real(8) p5_e,p5_l,p5_dl,p5_l2,p5_dl2
real(8) rxk1_w,rxk2_w
real(8) dr_well
real(8) Nideal_well(122),GRNideal_well(122),Gr_w(122)
real(8) hist_w(122),SumGr_w(122),VirGr_w(122),Grx_well(122)
real(8) sumP2,sumP2out,cv,cvout
real(8) kvec(10000)
integer(4) kmax,ksqmax
real(8) fmax
real(8) rmax(30,30,30,30)
!для модифицированных поетнциалов
real(8) kkm1(30,30,30,30),kkm2(30,30,30,30)
real(8) dkkm1(30,30,30,30),dkkm2(30,30,30,30)
real(8) modif
!
real(8) TorsCosFi,TorsFi
real(8) d0x,d0y,d0z
real(8) d1x,d1y,d1z
real(8) d2x,d2y,d2z
real(8) d0d1x,d0d1y,d0d1z
real(8) d1d2x,d1d2y,d1d2z
real(8),allocatable:: TorsDist(:,:,:)
!
integer(4) inter_tot
integer(4) inter_ac
integer(4) inter_nac
integer(4) intra_tot
integer(4) intra_nac
integer(4) intra_ac
real(8) RotDeg, TorsDeg, BindDeg
real(8) steptip
!
real(8) TotBindEn   !полная энергия напряженности валентных углов
real(8) TotTorsEn   !полная энергия торсионных углов
real(8) inter_tip, intra_tip
real(8) BindEn_do,BindEn_posle, deltaBind
real(8) TorsEn_do,TorsEn_posle, deltaTors
real(8) SumBindout,sumBind,sumTorsout,sumTors
real(8) bind_en_s,bind_en_out,tors_en_s,tors_en_out
real(8) sbind,stors
real(8) bind_out_out,tors_out_out
real(8) SumCvOut, OutCvOut
real(8) lp,Dep
real(8) DDLJatom
real(8) DDmol, DDVir,TotalDDLJ
real(8) DDLJ_do,DDLJ_posle,deltaDDLJ,sumDDLJ
real(8) sumVir2Out, betminus1,bette
real(8) o2cv
!---для трехчастичного взаимодействия
real(8) tr_cut      !обрезание трехчастичного потенциала
real(8),allocatable:: tr_rast(:,:)     !растояние от центра одной молекулы до центра другой
!
real(8),allocatable:: tri_x(:,:)   !Растояние до атомов без минимального образа
real(8),allocatable:: tri_y(:,:)
real(8),allocatable:: tri_z(:,:)
real(8),allocatable:: trx_l(:)
real(8),allocatable:: try_l(:)
real(8),allocatable:: trz_l(:)
real(8) Triple_en, Tr_en_do,Tr_en_posle
real(8) tr_v(10,3,10,3,10,3)
real(8) TotalTrEn
real(8) Uijk
integer(4) calc_tr
!модифициораванное перемещение
integer(4), allocatable:: near_mol(:)
real(8), allocatable:: near_rast(:)
real(8) min_rast_m
real(8) sumTrEn,tr_out,tr_dout
character(20) block_out
!
real(8) DbEn, DbVir, TotDbEn
integer(4) bl_no
real(8) glomin
integer(4) cnkol,ndismol
real(8) zgran1,zgran2, zdel
real(8), allocatable:: denssum(:),denssum1(:),denssum2(:)
integer(4) densl,densl1,densl2
integer(4) densNo
end module dannie
!****************************************************************************
module generator
 integer(4) n1, n2, n3 !случаные числа
 integer(4) m1,a1,b1,a2,b2,m3,m2 !параметры генератора
 integer(4) max1, max2,max3 !максимумы генераторов
 real(8) randmass(500) !массив случайных чисел
 real(8) outrand
end module
!Изменение длинны шага

subroutine chDlShag !динамическое изменение длинны шага
use dannie

if (intra_tot>0) then      !изменение максимального шага
  !print *, intra_ac, intra_tot, float(intra_ac)/float(intra_tot)
 if (float(intra_ac)/float(intra_tot)<0.4) then
  DlShag=DlShag/1.1
  RotDeg=RotDeg/1.1
 endif
 if (float(intra_ac)/float(intra_tot)>0.6) then
  DlShag=DlShag*1.1
  RotDeg=RotDeg*1.1
 endif
endif
if (DlShag>KonSt*0.8) then       !ограничение диапазона изменения шага
 DlShag=KonSt*0.8
endif
if (DlShag<KonSt*0.005) then
 DlShag=KonSt*0.005
endif
if (RotDeg>0.2) then          !максимальный угол поворота по одной оси 17 градусов
 RotDeg=0.2
endif
if (RotDeg<0.1) then          !минимальный угол поворота по одной оси 5 градусов
 RotDeg=0.1
endif
RotDeg=0.5
if (inter_tot>0) then      !изменение максимального изменения
 if (float(inter_ac)/float(inter_tot)<0.4) then
  TorsDeg=TorsDeg/1.1
 endif
 if (float(inter_ac)/float(inter_tot)>0.6) then
  TorsDeg=TorsDeg*1.1
 endif
endif
if (TorsDeg<0.07) then   !ограничение максимального изменения угла в диапазоне
 TorsDeg=0.07            !от 4 до 30 градусов
endif
if (TorsDeg>0.5) then
 TorsDeg=0.5
endif
!
end subroutine chDlShag
!****************************************************************************
subroutine input
use dannie
real(8) ix(100)
real(8) iy(100)
real(8) iz(100)
real(8) ir(100)
real(8) iep(100)
real(8) isi(100)
real(8) ia(100)
real(8) ixc(100)
real(8) iyc(100)
real(8) izc(100)
character(20) iname
real(8) MMi
integer(4) Nmi
character(5) ill(100)
character(5) ilc(100)
integer(4) torsi !количество торсионных связей в молекуле
integer(4) fsti(30), seci(30), nroti(100)
integer(4) rotai(100,100)
integer(4) srotai(100,100)
real(8) torskoefi(100,5)
integer(4) NMCI
real(8) torsCosi(30)
integer(4) bind_koli
integer(4) bind_ci(30)
integer(4) bind1i(30)
integer(4) bind2i(30)
integer(4) bind1ki(30,30)
integer(4) bind2ki(30,30)
real(8) bind_ki(30)
real(8) bind_fii(30)
integer(4) l_koli(30)
real(8) l_ki(30)
real(8) l_ri(30)
integer(4) l_nk1i(30,30)
integer(4) l_nk2i(30,30)
integer(4) bind_vali
integer(4) bind_vni(30,3)
integer(4) bond_koli
integer(4) bond_ni(30,2)
integer(4) charg_koli
real(8) ch_xi(30)
real(8) ch_yi(30)
real(8) ch_zi(30)
real(8) ch_qi(30)
integer(4) ch_bind1i(30)
integer(4) ch_bind2i(30)
integer(4) ch_bind1ki(30,30)
integer(4) ch_bind2ki(30,30)
integer(4) ti,tia,tj,tja,tz,tza
!****************************************************************************
!Ввод исходных данных
namelist /main/RoNV,Nv,sproc,storona,Naz,tipResh,Temp,RadCut,TipCut, TipLat&
&,Nchan,drtail,Nstep,drtail,workdir,ena,npr,ptip,calc_rdf,calc_hend, calc_tr
!нужно использовать временные координаты
namelist /molekula/MMi,Nmi,iName,ix,iy,iz,ir,iep,ill,ia,isi,nmci,ixc,iyc&
&,izc,ilc,torsi, fsti, seci, nroti, rotai,srotai, torskoefi, &
& bind_koli,bind_ci, bind1i,bind1ki,bind2i,bind2ki,bind_ki&
&,bind_fii, l_koli,l_ki,l_ri,l_nk1i,l_nk2i, bind_vali,bind_vni,&
&bond_koli,bond_ni, charg_koli,ch_xi,ch_yi,ch_zi,ch_qi,ch_bind1i&
&,ch_bind2i,ch_bind1ki,ch_bind2ki !данные молекулы (параметры атомов)
!namelist /config/Neqv,Ntot,Nnew,Nchan !всего перемещений, перемещений в записи, смена шага
open(19,file='input.txt')
 read(19,main)
close(19)
!
!lentr=len_trim(workdir)
workdir=adjustr(workdir)
!print *, workdir
do i=1,Nv
 Nazv=Naz(i)
 write(Nazv2,'(a4,a30)') '/VV/',Nazv
 write(NazF(i),'(a200,a24)') workdir, Nazv2
 NazF(i)=adjustl(NazF(i))
enddo
do i=1,Nv,1
 open(10,file=NazF(i))
  read(10,molekula) !считывание данных по молекуле
  MM(i)=MMi
  Nm(i)=Nmi
  Mcm(i)=NMCI !количество атомов
  MMName(i)=iName
  do j=1,100,1
   xm(i,j)=ix(j)
   ym(i,j)=iy(j)
   zm(i,j)=iz(j)
   rm(i,j)=ir(j)
   xcm(i,j)=ixc(j) !координаты атомов
   ycm(i,j)=iyc(j)
   zcm(i,j)=izc(j)
   labm(i,j)=ilc(j)
   epsi(i,j)=iep(j)
   sigma(i,j)=isi(j)
   alfa(i,j)=ia(j)
   labelm(i,j)=ill(j)
  enddo
  !********Торсионный угол******
  Ntors(i)=torsi
   do j=1,Ntors(i),1!если есть торсионнные связи
    fsta(i,j)=fsti(j) !номер нулевого атома при повороте
    seca(i,j)=seci(j) !номер атома на оси
    nrota(i,j)=nroti(j) !количество поворачиваемых частиц
    do k=1,nroti(j),1 !цикл считывания поворочиваемых центров
     rota(i,j,k)=rotai(j,k) !номера поворачиваемых частиц
     srota(i,j,k)=srotai(j,k) !номера частиц около второго атома
    enddo
    do k=1,5,1
     torskoef(i,j,k)=torskoefi(j,k)
    enddo
   enddo
   !*******Валентный угол******
   Bind_kol(i)=Bind_koli !количество валентных углов
   do j=1,bind_kol(i),1
    bind_cent(i,j)=bind_ci(j)        !номер центра
    bind1(i,j)=bind1i(j)             !количество центров с одной стороны
    bind2(i,j)=bind2i(j)             !количество центров с другой
    ch_bind1(i,j)=ch_bind1i(j)
    ch_bind2(i,j)=ch_bind2i(j)
    !print *,bind_cent(i,j),bind1(i,j),bind2(i,j)
    do k=1,bind1(i,j),1
     bind1k(i,j,k)=bind1ki(j,k)      !номера поворачиваемых центров
     !первый центр соотвествует центру составляющему угол
     !print *,bind1k(i,j,k)
    enddo
    do k=1,bind2(i,j),1
     bind2k(i,j,k)=bind2ki(j,k)      !номера поворачиваемых центров
     !первый центр соотвествует центру составляющему угол
     !print *, bind2k(i,j,k)
    enddo
    do k=1,ch_bind1(i,j),1
     ch_bind1k(i,j,k)=ch_bind1ki(j,k)
     ch_bind2k(i,j,k)=ch_bind2ki(j,k)
    enddo
    bind_k(i,j)=bind_ki(j)/2.0
    bind_fi(i,j)=bind_fii(j)/180.0*ChPI
    !*******Длинна связи********
    !l_kol
    !do j=1,
    !print *,bind_fi(i,j),bind_k(i,j)
   enddo
   bind_val(i)=bind_vali         !валентный угол
   !print *, bind_vali
   !pause
   do j=1,bind_val(i)
    bind_vn(i,j,1)=bind_vni(j,1)
    bind_vn(i,j,2)=bind_vni(j,2)
    bind_vn(i,j,3)=bind_vni(j,3)
   enddo
   bond_kol(i)=bond_koli
   do j=1,bond_kol(i),1
    bond_n(i,j,1)=bond_ni(j,1)
    bond_n(i,j,2)=bond_ni(j,2)
   enddo
   !заряды
 !РАсставляем заряды на молекулах
  charg_kol(i)=charg_koli
  do j=1,charg_kol(i),1
   ch_x(i,j)=ch_xi(j)
   ch_y(i,j)=ch_yi(j)
   ch_z(i,j)=ch_zi(j)
   ch_q(i,j)=ch_qi(j)
  enddo
  close(10)
enddo
! определяем количество расставляемых молекул
select case (tipResh)
case (2)
 N=4*(storona)**3 !определяем количество молекул
case (1)
 N=2*(storona)**3
end select
!
!определяем безразмерный коэффициент
bs=1.0
be=1.0
do i=1,Nv,1
 do j=1,Nm(i),1
  sigma(i,j)=sigma(i,j)/bs
  epsi(i,j)=epsi(i,j)/be
  xm(i,j)=xm(i,j)/bs
  ym(i,j)=ym(i,j)/bs
  zm(i,j)=zm(i,j)/bs
  if (tr_cut<sigma(i,j)) then
   tr_cut=sigma(i,j)
  endif
  !print *, xm(i,j),ym(i,j),zm(i,j)
 enddo
enddo
tr_cut=tr_cut*2.0
!pause
KoefP=0.3 !коэффициент поворота
if (ena/=1) then
 bezro=bs*bs*bs*6.022045/10000.0
 bezt=1/be
 beze=be*8.3145107
 bezp=1/(bs*bs*bs)*be*100.0*1.38065812
 bezcv=8.3145107
 razm_davl=' bar '
 razm_ro=' mol/liter '
 razm_tem=' K '
 razm_en=' J/mol'
 razm_rast=' Å '
 razm_mu=' J/mol '
 razm_cv=' J/(mol·K) '
else
 bezro=bs*bs*bs
 bezt=1/be
 beze=be
 bezp=1/(bs*bs*bs)*be
 bezcv=1
 razm_davl=' '
 razm_ro=' '
 razm_tem=' '
 razm_en=' '
 razm_rast=' '
 razm_mu=' '
 razm_cv=' '
endif
RoNVML=RoNV
TempK=Temp
RoNV=RoNV*bezro!*bs*bs*bs*6.022045/10000.0
Temp=Temp*bezt!/be
!Закончили ввод данных
!Параметры трехчастичного взаимодействия
do ti=1,Nv,1
 do tia=1,Nm(ti),1
  do tj=1,Nv,1
   do tja=1,Nm(tj),1
    do tz=1,Nv,1
     do tza=1,Nm(tz),1
      tr_v(ti,tia,tj,tja,tz,tza)=3.0
     enddo
    enddo
   enddo
  enddo
 enddo
enddo

end subroutine

!****************************************************************************
program mcprog
!INCLUDE "omp_lib.h"
use dannie
integer(4) Nmol !номер передвигаемой/поворачиваемой молекулы
integer(4) i,j,k,l,mi
!Закончили объявление переменных
!****************************************************************************
call cpu_time(start)
ChPI=3.14159265358979323846264338327
ChPI2=2.0*ChPI
!ДЛЯ суммирования эвалда
MRDF=0
MRDFN=0
kmax=5
ksqmax=27
modif=7            !ОПРЕДЕЛЯТЬ ЗДЕСЬ
glomin=0.6
cnkol=0
ndismol=0
open(59,file='movie.xyz')
open(56,file='ndis.txt')
!!$OMP PARALLEL
!!k=OMP_GET_NUM_PROCS()
!!print *, k
!!$OMP END PARALLEL
!!pause
allocate(kolNv(Nv+1))
versia='1.02.8'
call input !ввод исходных данных

call randomn !задействование случайных чисел
!do i=1,50,1
! print *, getrand()
!enddo
call newlatice !расстановка молекул в объем (получение параметрова рещетки)
!проверка на нужность вращения и тд
inter_tip=1.0
intra_tip=0.0
do i=1,Nv,1
 if (Ntors(i)>0) then
  inter_tip=0.5
  intra_tip=0.0
 endif
enddo
!write(6,'(f4.2,a,f4.2,a,f4.2)') torsch_tip,'  -  ', rotch_tip,'  -  '&
!&,rotch_tip
RotDeg=0.25
TorsDeg=0.09
BindDeg=0.02
if (calc_hend==1) then
 NChan=N*32*2*2*2    !*2 !*2
 NNew=NChan*8 !ло просто *16  !2048*N !2048*N   !было 2048
 Neqv=4*Nnew     !было 3
 Ntot=Neqv
else
 NChan=N*32*2*2*2    !*2 !*2
 if (calc_tr==1) then
  NChan=NChan
 endif
 NNew=NChan*16 !ло просто *16  !2048*N !2048*N   !было 2048
 Neqv=4*Nnew     !было 3
 Ntot=4*Neqv
 !!новые значения с большым пробегом на равновесие
 !Neqv=30*1000000
 !NChan=100000
 !NNew=10000000
 !Ntot=NNew*8 !!
endif
!!считывание конфигурации из файла
allocate(oenp(Ntot/Nnew))
allocate(odavl(Ntot/Nnew))
allocate(ocv(Ntot/Nnew))
allocate(xaado(NmMax+1)) !абсолютные координаты
allocate(yaado(NmMax+1)) !
allocate(zaado(NmMax+1)) !
allocate(xado(NmMax+1)) !абсолютные координаты
allocate(yado(NmMax+1)) !
allocate(zado(NmMax+1)) !
allocate(TorsDist(Nv,NtMax,720))
!ЧТЕНИЕ ИЗ ФАЙЛА !!!ПРОДУМАТЬ
!if (TipLat==2) then
!open (12,file='INLATM.txt')
!do i=1,N,1
! read (13,'(3f30.20)') x(i),y(i),z(i)
! x(i)=x(i)*konst
! y(i)=y(i)*konst
! z(i)=z(i)*konst
!enddo
!close (12)
!open (13,file='INLATA.txt')
! do i=1,Na,1
!  read (13,'(3f30.20)') xa(i),ya(i),za(i)
! enddo
!close (13)
!endif
!Заводим матрицы для сигмы
!allocate(siga(Nv,Nv,NmMax,NmMax))
!allocate(epsa(Nv,Nv,NmMax,NmMax))
!обуление файла шагов
 open(17,file='OUTstep.txt')
 write (17,'(a)') 'steps; energy; average energy; standard &
 & deviation; pressure; average pressure; standard deviation'
 close(17)

do i=1,Nv,1
 do j=1,Nv,1
  do k=1,Nm(i),1
   do l=1,Nm(j),1
   !для леннарда джонса
    epsa(i,j,k,l)=sqrt(epsi(i,k)*epsi(j,l))
    alf(i,j,k,l)=sqrt(alfa(i,k)*alfa(j,l))
    siga(i,j,k,l)=(sigma(i,k)+sigma(j,l))/2.0
    if (ptip==9) then
     dopeps(i,j,k,l)=(alf(i,j,k,l)/(alf(i,j,k,l)-6))*(alf(i,j,k,l)/6.0)**(6.0/(alf(i,j,k,l)-6.0))
     !print *, dopeps(i,j,k,l)
    endif
    siga6=siga(i,j,k,l)**6
    !карра-коновалова
    kk1(i,j,k,l)=6.0*epsa(i,j,k,l)*siga6/alf(i,j,k,l)
    kk2(i,j,k,l)=-epsa(i,j,k,l)*(alf(i,j,k,l)+6.0)*siga6/alf(i,j,k,l)
    !
    kkm1(i,j,k,l)=modif*epsa(i,j,k,l)*siga(i,j,k,l)**modif/alf(i,j,k,l)
    kkm2(i,j,k,l)=-epsa(i,j,k,l)*(alf(i,j,k,l)+modif)*siga(i,j,k,l)**modif/alf(i,j,k,l)
    !print *,kkm1(i,j,k,l),kkm2(i,j,k,l)
    !pause
    !kk3(i,j,k,l)=-6.0*epsa(i,j,k,l)*siga6/alf(i,j,k,l)
    !kk4(i,j,k,l)=-epsa(i,j,k,l)*siga6
    dkk1(i,j,k,l)=-6.0*epsa(i,j,k,l)*siga6/siga(i,j,k,l)
    dkk2(i,j,k,l)=-36.0*epsa(i,j,k,l)*siga6/alf(i,j,k,l)
    dkk3(i,j,k,l)=-dkk2(i,j,k,l)
    dkk4(i,j,k,l)=6.0*epsa(i,j,k,l)*siga6
    dkkm1(i,j,k,l)=-epsa(i,j,k,l)*modif*(modif+alf(i,j,k,l))/alf(i,j,k,l)*siga(i,j,k,l)**modif
    dkkm2(i,j,k,l)=-epsa(i,j,k,l)*modif/siga(i,j,k,l)*siga(i,j,k,l)**modif
    !print *, kk1(i,j,k,l),kk2(i,j,k,l),kk3(i,j,k,l),kk4(i,j,k,l)
    !для букингема
    bk1(i,j,k,l)=6.0*epsa(i,j,k,l)/(alf(i,j,k,l)-6.0)
    bk2(i,j,k,l)=-epsa(i,j,k,l)*siga6/(1.0-6.0/alf(i,j,k,l))
    dbk1(i,j,k,l)=6.0*epsa(i,j,k,l)*siga6/(1.0-6.0/alf(i,j,k,l))
    dbk2(i,j,k,l)=-6.0*epsa(i,j,k,l)/(siga(i,j,k,l)-6.0*siga(i,j,k,l)/alf(i,j,k,l))
    fmax=-100000
    if (ptip==4) then
     do zi=1,1000,1
      rz=zi*siga(i,j,k,l)/1000.0
      rz2=1.0/rz/rz
      rz6=rz2*rz2*rz2
      expp=exp(alf(i,j,k,l)*(1.0-rz/siga(i,j,k,l)))
      ekt=bk1(i,j,k,l)*expp+bk2(i,j,k,l)*rz6
      if (ekt>fmax) then
       rmax(i,j,k,l)=rz
       fmax=ekt
      endif
     enddo
     write(6,'(a, f15.10)') 'Buckingham maximum:  ', rmax(i,j,k,l)
     !rmax(i,j,k,l)=rmax(i,j,k,l)  !*rmax(i,j,k,l)
    endif
    !print *, bk1(i,j,k,l),bk2(i,j,k,l)
   enddo
  enddo
 enddo
enddo
!Если потенциал Букингема -определяем максимальное растояние
call creatAA !растановка атомов в объем
!вывод исходных данных
write (6,'(a80)') '**********************************MCPROG***************************************'
write (6,'(a,a)') 'Program internal version: ', trim(adjustl(versia))
if (ptip==1) then
 write (6,'(a)') trim(adjustl('Potential:   Karr-Konowalow'))
endif
if (ptip==2) then
 write (6,'(a)') trim(adjustl('Potential:   Lennard-Jones (12-6)'))
endif
if (ptip==3) then
 lp=0.3     !установка параметров
 Dep=0.0
 write (6,'(a)') trim(adjustl('Potential:   Lennard-Jones+Potential Well'))
 write (6,'(a,f15.5)') trim(adjustl('l*: ')),lp
 write (6,'(a,f15.5)') trim(adjustl('De: ')),Dep
endif
if (ptip==4) then
 write (6,'(a)') trim(adjustl('Potential:   Buckingham (exp-6)'))
endif
if (ptip==5) then
 p5_e=0.0                  !установка пармтеров
 p5_l=0.3
 p5_dl=0.025
 p5_l2=p5_l*p5_l
 p5_dl2=(p5_l+p5_dl)*(p5_l+p5_dl)
 write (6,'(a)') 'Potential: Karr-Konowalow + square-well'
 write (6,'(a,f15.7)') 'e*: ', p5_e
 write (6,'(a,f15.7)') 'l:  ', p5_l
 write (6,'(a,f15.7)') 'l+Δl', p5_l+p5_dl
endif
if (ptip==7) then
 write(6,'(a)') 'Potential: Modifed Karr-Konowalow'
 write(6,'(a,f6.4)') 'Attractive parameter: ', modif
endif
if (ptip==9) then
 write(6,'(a)') 'Potential: n-6 Lennard-Jones'
endif
write(6,'(a)') '______________Molecular properties______________'
do i=1,Nv,1
 write (6,'(a,a,a3,i3)') 'Centers in ',trim(adjustl(MMName(i))),': ' ,Nm(i)
 write (6,'(a,a,a3,f7.3,a,i4,a)') 'Mol. fraction of ',trim(adjustl(MMName(i))),': ',&
 & Sproc(i)*100.0,' % ',kolNv(i), ' particles'
 if ((ptip==1).or.(ptip==7)) then !потенциал Карра -коновалова
  do j=1,Nm(Nv),1
   write (6,'(a,a,f6.3,a,f7.4,a,a,f10.5,a)') trim(adjustl(labelm(i,j))), ' α ', alfa(i,j),&
   &' rm ', sigma(i,j), trim(razm_rast), &
   & ' ε/kb ', epsi(i,j), trim(razm_tem)
  enddo
 endif
 if ((ptip==2).or.(ptip==5) .or. (ptip==3)) then !потенциал Леннард-Джонса
  do j=1,Nm(Nv),1
   write (6,'(a,a,f7.4,a,a,f10.5,a)') trim(adjustl(labelm(i,j))),' σ ', sigma(i,j),&
   & trim(razm_rast),' ε/kb ', epsi(i,j), trim(razm_tem)
  enddo
 endif
  if (ptip==4) then !потенциал Букингема
  do j=1,Nm(Nv),1
   write (6,'(a,a,f6.3,a,f7.4,a,a,f10.5,a)') trim(adjustl(labelm(i,j))), ' α ', alfa(i,j),&
   &' rm ', sigma(i,j), trim(razm_rast), &
   & ' ε/kb ', epsi(i,j), trim(razm_tem)
  enddo
 endif
 if (ptip==9) then
  do j=1,Nm(Nv),1
   write (6,'(a,a,f6.3,a,f7.4,a,a,f10.5,a)') trim(adjustl(labelm(i,j))), ' n ', alfa(i,j),&
   &' σ ', sigma(i,j), trim(razm_rast), &
   & ' ε/kb ', epsi(i,j), trim(razm_tem)
  enddo
 endif
 do j=1,ckol,1
  write (6,'(a,a3,i5)') trim(adjustl(labckol(j))),' - ', Ncrdf(i,j)
 enddo
 !длинны связи
 if (bond_kol(i)>0) then
  write(6,'(a)') 'Start bond length: '
  do j=1, bond_kol(i),1
   x1valent=xm(i,bond_n(i,j,1))-xm(i,bond_n(i,j,2))
   y1valent=ym(i,bond_n(i,j,1))-ym(i,bond_n(i,j,2))
   z1valent=zm(i,bond_n(i,j,1))-zm(i,bond_n(i,j,2))
   write(6,'(a,a,a,a,f10.5,a)') trim(adjustl(labelm(i,bond_n(i,j,1)))), '-',trim(adjustl(labelm&
   &(i,bond_n(i,j,2)))), '  : ', sqrt(x1valent*x1valent+&
   &y1valent*y1valent+z1valent*z1valent), razm_rast
  enddo
 else
  write (6,'(a)') 'Chemical bonds: None '
 endif
 !валентные углы
 if (bind_val(i)>0) then
 write(6,'(a)') 'Start valent angles: '
 do j=1,bind_val(i),1
  x1valent=xm(i,bind_vn(i,j,1))-xm(i,bind_vn(i,j,2))
  x2valent=xm(i,bind_vn(i,j,3))-xm(i,bind_vn(i,j,2))
  y1valent=ym(i,bind_vn(i,j,1))-ym(i,bind_vn(i,j,2))
  y2valent=ym(i,bind_vn(i,j,3))-ym(i,bind_vn(i,j,2))
  z1valent=zm(i,bind_vn(i,j,1))-zm(i,bind_vn(i,j,2))
  z2valent=zm(i,bind_vn(i,j,3))-zm(i,bind_vn(i,j,2))
  valcos=(x1valent*x2valent+y1valent*y2valent+z1valent*z2valent)/&
  &sqrt(x1valent*x1valent+y1valent*y1valent+z1valent*z1valent)/&
  &sqrt(x2valent*x2valent+y2valent*y2valent+z2valent*z2valent)
  !print *, sqrt(x1valent*x1valent+y1valent*y1valent+z1valent*z1valent)
  !print *, sqrt(x2valent*x2valent+y2valent*y2valent+z2valent*z2valent)
  !расчет валентных углов
  write(6,'(6a,f12.7,a)') trim(adjustl(labelm(i,bind_vn(i,j,1)))), '-',&
  trim(adjustl(labelm(i,bind_vn(i,j,2)))), '-',  &
  &trim(adjustl(labelm(i,bind_vn(i,j,3)))), ' :',  acos(valcos)/ChPI*180.0, ' º'
 enddo
 else
  write (6,'(a)') 'Valent angles: None '
 endif
 if (bind_kol(i)>0) then
  write(6, '(a)') 'Valent angle not rigid'
  do j=1,bind_kol(i),1
   write(6,'(6a)') trim(adjustl(labelm(i,bind1(i,j)))), '-',trim(adjustl(labelm(&
   &i,bind_cent(i,j)))),'-',trim(adjustl(labelm(i,bind2(i,j)))), ' :'
   write(6, '(a,f20.10)') '                       k : ', bind_k(i,j)*2.0
   write(6, '(a,f20.10)') '       eqvilibrium angle : ', bind_fi(i,j)/ChPI*180.0
  enddo
 endif
 !Торсионные углы
 write(6,'(a)') 'Start torsion angles:'
 do j=1,Ntors(i),1
  write(6,'(4a)') trim(adjustl(labelm(i,fsta(i,j)))),'-',trim(adjustl(labelm(i,seca(i,j))))&
  &,'  : '
  write(6,'(a,i2,a)') 'Rotating sites (', Nrota(i,j), ')'
  do k=1,Nrota(i,j)
   write(6,'(a,a,i2)')'     ', trim(adjustl(labelm(i,rota(i,j,k)))),rota(i,j,k)
  enddo
  do k=1,3,1
   write(6,'(a,i1,a,f20.10)') 'V',k,' :  ', TorsKoef(i,j,k)
  enddo
  !
   !определение торсионого угла
   d0x=xm(i,rota(i,j,1))-xm(i,fsta(i,j))
   d0y=ym(i,rota(i,j,1))-ym(i,fsta(i,j))
   d0z=zm(i,rota(i,j,1))-zm(i,fsta(i,j))
   !
   d1x=xm(i,seca(i,j))-xa(i,fsta(i,j))
   d1y=ym(i,seca(i,j))-ya(i,fsta(i,j))
   d1z=zm(i,seca(i,j))-za(i,fsta(i,j))
   !
   d2x=xa(i,srota(i,j,1))-xa(i,seca(i,j))
   d2y=ya(i,srota(i,j,1))-ya(i,seca(i,j))
   d2z=za(i,srota(i,j,1))-za(i,seca(i,j))
   !
   d0d1x=d0y*d1z-d0z*d1y
   d0d1y=d0z*d1x-d0x*d1z
   d0d1z=d0x*d1y-d0y*d1x
   !
   d1d2x=d1y*d2z-d1z*d2y
   d1d2y=d1z*d2x-d1x*d2z
   d1d2z=d1x*d2y-d1y*d2x
   !
   TorsCosFi=(d0d1x*d1d2x+d0d1y*d1d2y+d0d1z*d1d2z)/sqrt(d0d1x*d0d1x+d0d1y*d0d1y+&
   &d0d1z*d0d1z)/sqrt(d1d2x*d1d2x+d1d2y*d1d2y+d1d2z*d1d2z)
   TorsFi=acos(TorsCosFi)/ChPI*180.0
   if (TorsCosFi==1) then
    if (d0d1x*d1d2x+d0d1y*d1d2y+d0d1z*d1d2z>0.00000001) then
     TorsFi=180.0
    endif
   endif
   write(6,'(a,f20.10,a,f20.10,a)') 'start angle cos: ', TorsCosFi, 'angle:  ',  &
   &TorsFi, ' º'
   !pause
  !
 enddo
 !
  !ЗАряды

 if (charg_kol(i)>0) then
  write(6,'(a)') 'Charges :'
  if (charg_kol(i)>0) then
   do j=1,charg_kol(i),1
    write(6,'(f15.10)') ch_q(i,j)
   enddo
  endif
 else
  write (6,'(a)') 'Charges: None'
 endif
 !Длины связи


 write (6,'(a)') '_______________________________________________________'
enddo
write (6,'(a,i5)') trim(adjustl('Molecules in volume: ')), N
write (6,'(a,f10.4,a)') trim(adjustl('Mol. density: ')),RoNVML, trim(razm_ro)
write (6,'(a,f10.4,a2)') trim(adjustl('Temperature: ')), TempK, trim(razm_tem)
write (6,'(a,i10)') trim(adjustl('Equilibrium steps: ')), Neqv
write (6,'(a,i10)') trim(adjustl('Sample steps: ')), Nnew
write (6,'(a,i10)') trim(adjustl('Total steps: ')), Ntot
!if (TipCut==1) then
! write (6,'(a30)') 'Tip obrezania mol-mol'
!else
! write (6,'(a30)') 'Tip obrezania atom-atom'
!end if
if (TipLat==1) then
 if (tipResh==1) then
 write (6,'(a)') trim(adjustl('Start lattice - OCK'))
 else if (tipResh==2) then
 write (6,'(a)') trim(adjustl('Start latice - GCK'))
 end if
end if
do i=1,ckol,1
 write (6,'(i3,a,a)') i, ' - ', labckol(i)
enddo
write (6,'(a)')'________________________________________________________'
!Проверка соразмерности ячейки и радиуса обрезания
if (RadCut<KonSt/2) then
 write (6,'(a,f15.8,a)') trim(adjustl('Cut radius: ')), RadCut, trim(razm_rast)
else
 RadCut=KonSt/2
 write (6,'(a,f15.8,a)') trim(adjustl('Cut radius: ')), RadCut, trim(razm_rast)
endif
call tails() !Определяем хвостовой потенциал
RadCut2=RadCut*RadCut !квадрат растояния обрезания
rGrcut=KonSt/2.0 !поставить половину длинны объема
stepen=8
deltar=100.0
!
!open(58,file='starttors.txt')
!do i=1,N,1
! do j=1,Ntors(Ntip(i))
!  call calc_tors(i,j)
!  write (58,'(i5,i5,f20.10,f20.10)') i,j,TorsCosFi,TorsFi
! enddo
!enddo
!close(58)
!!pause
!!
!максимальная сигма
do i=1,Nv,1
 if (maxsig<sigma(i,1)) then
  maxsig=sigma(i,1)
 endif
enddo
do while (deltar>0.06*maxsig)     !раньше было 0.005
 stepen=stepen+1
 HistRaz=2**stepen !2**11 !на один меньше (добавим нулевой в конце)
 deltar=(rGrcut)/(float(HistRaz)+0.5)
enddo
!
!kol_m=ceiling(20.0/RadCut)-1
kol_m=2
rGrcut_m=rGrcut*kol_m
N_m=N*kol_m*kol_m*kol_m
HistRaz_m=kol_m*HistRaz
Mrdf_m=0
allocate(x_m(N_m))
allocate(y_m(N_m))
allocate(z_m(N_m))
allocate(xaa_m(Nv,N_m))
allocate(yaa_m(Nv,N_m))
allocate(zaa_m(Nv,N_m))
allocate(HistAtom_m(kolfrr,HistRaz_m))
allocate(GrA_m(kolfrr,HistRaz_m))
allocate(VirGrA_m(kolfrr,HistRaz_m))
allocate(smooth(Histraz+10))
allocate(Grx_m(HistRaz_m))
allocate(Nideal_m(HistRaz_m))
allocate(total_c_f_m(kolfrr,HistRaz_m))
allocate(SumGrA_m(kolfrr,HistRaz_m))
allocate(Bridge_f(kolfrr,Histraz))
!
!write (6,'(a,i6,a,f15.8,a)') 'cell number ', HistRaz,  ' cell length ', deltar,trim(razm_rast)
!write (6,'(a,i6)') 'addition volume number', kol_m
delta_rast=deltar
kon_rast=HistRaz
allocate(HistY(Nv,HistRaz))
allocate(Hist(HistRaz)) !полос гистограммы
allocate(Gr(HistRaz)) !координат распределения
allocate(Nideal(HistRaz)) !распределение для идеального газа
allocate(Grx(HistRaz))
allocate(SumGr(HistRaz))
allocate(VirGr(HistRaz))
allocate(HistAtom(kolfrr,HistRaz))
allocate(GrA(kolfrr,HistRaz))
allocate(SumGrA(kolfrr,HistRaz))
allocate(VirGrA(kolfrr,HistRaz))
allocate(SumGrAN(kolfrr,HistRaz))
allocate(VirGrAN(kolfrr,HistRaz))
allocate(RDF_Block(20,kolfrr,HistRaz))
allocate(total_c_f(kolfrr,HistRaz))
allocate(direct_c_f(kolfrr,HistRaz))
allocate(indirect_c_f(kolfrr,HistRaz))
allocate(sum_minus(Nv,HistRaz))
allocate(y_potenc2(Nv,kon_rast))
allocate(cavity_c_f(Nv,kon_rast))
allocate(sum_potenc2(Nv,kon_rast))
allocate(cavity_out(Nv,HistRaz))
allocate(nom_y(Nv,HistRaz))
allocate(mu2_nom(Nv))
allocate(mu2_sum(Nv))
allocate(mu2_out(Nv))
! for FFT
allocate(AFFT(HistRaz))
allocate(BFFT(Histraz))
allocate(WorkA(HistRaz*2))
allocate(WorkB(HistRaz*2))
!
do i=1,HistRaz,1
 SumGr(i)=0.0
 VirGr(i)=0.0
enddo
do inom=1,ckol,1
 do jnom=inom,ckol,1
  do k=1,HistRaz,1
   SumGrA((inom-1)*ckol+jnom,k)=0.0
   VirGrA((inom-1)*ckol+jnom,k)=0.0
  enddo
 enddo
enddo
GRconst=4.0*ChPi*RoNV/3.0
!находим количество атомов в молекуле
atomk=0
do i=1,N,1
 do j=1,Na(i),1
  atomk=atomk+1
 enddo
enddo
!нахождение идеального распределения и растояния
do i=1,(HistRaz),1
 rNiz=(i-0.5)*deltar
 rVerh=(i+0.5)*deltar
 Nideal(i)=GRconst*(rVerh*rVerh*rVerh-rNiz*rNiz*rNiz)
 Grx(i)=float(i)*deltar
 !if ((Grx(i)>0.3).and.(Grx(i)<0.35)) then
 ! print *,Grx(i),Nideal(i)
 !endif
enddo
!нахождение идеального распеределения для потенциалов с ямой

if (ptip==5) then            !прямоугольная потенциальная яма
 dr_well=p5_dl/120.0
 do i=1,122,1
  rNiz=p5_l+float(i-2)*dr_well
  rVerh=p5_l+float(i-1)*dr_well
  Nideal_well(i)=GRconst*(rVerh*rVerh*rVerh-rNiz*rNiz*rNiz)
  Grx_well(i)=(rNiz+rVerh)/2.0
  !print *, Grx_well(i),Nideal_well(i)
 enddo
 rxk1_w=p5_l-dr_well                !границы искомого фрр
 rxk2_w=p5_l+p5_dl+dr_well
 !print *, rxk1_w,rxk2_w,dr_well
 !pause
endif
if (ptip==3) then
 dr_well=0.002
 do i=1,122,1
  rNiz=lp+float(i-62)*dr_well
  rVerh=lp+float(i-61)*dr_well
  Nideal_well(i)=Grconst*(rVerh*rVerh*rVerh-rNiz*rNiz*rNiz)
  Grx_well(i)=(rNiz+rVerh)/2.0
  !print *, Nideal_well(i),Grx_well(i)
 enddo
 rxk1_w=lp-dr_well*61.0           !границы искомого ФРР
 rxk2_w=lp+dr_well*61.0
 !print *, rxk1_w,rxk2_w
 !pause
endif
!нахождение кинетической энергии
do i=1,(HistRaz_m),1
 rNiz=(i-1.0)*deltar
 rVerh=i*deltar
 Nideal_m(i)=GRconst*(rVerh*rVerh*rVerh-rNiz*rNiz*rNiz)
 Grx_m(i)=(i-0.5)*deltar
enddo
enid=0.0
do mi=1,Nv,1
 mu2_nom(mi)=0 !бнуление для нахождения химического потенциала
 mu2_out(mi)=0.0 !
 mu2_sum(mi)=0.0
 if (Nm(mi)==1) then
  enid=3.0/2.0*temp*sproc(mi)+enid
 else if (Nm(mi)==2) then
  enid=5.0/2.0*temp*sproc(mi)+enid
 else
  enid=3.0*temp*sproc(mi)+enid
 endif
enddo
! обнуление всего
assept=0
notassept=0
Nmow=0
Nout=0
naDl=0
aDl=0
sumP=0.0
sumVir=0.0
inter_tot=0
inter_ac=0
inter_nac=0
intra_tot=0
intra_nac=0
intra_ac=0
SumCvOut=0.0
OutCvOut=0.0
 sumPout=0.0
 sumP2out=0.0
 sumVirout=0.0
 sumVir2out=0.0
 sumBindout=0.0
 sumTorsout=0.0
 sumDDLJ=0.0
 sumTrEn=0.0
!call outrez    !проверка начала
!write (6,'(4(a,a6))') 'Step #                 E ',razm_en,'   P ', razm_davl,'   <E> ',razm_en,'   <P> ',razm_davl! вывод результатов на экран
!print *, '123123123123'
call totalen
!print *, '123123123123'
!print *, toten, totalvir
!   open(13,file='sumP.txt')
!проверка изменения валентного угла
zgran1=0.0
zgran2=KonSt
zdel=zgran2-zgran1
print *,'zdel',zdel
densl=90
densl1=60
densl2=30
allocate(denssum(densl))
allocate(denssum1(densl1))
allocate(denssum2(densl2))
do i=1,densl
    denssum(i)=0.0
enddo
do i=1,densl1
    denssum1(i)=0.0
enddo
do i=1,densl2
    denssum2(i)=0.0
enddo
do while (mov<=Ntot+Neqv) !начало шага монте карло
 mov=mov+1 !плюс шаг
 randn=getrand() !выбор случайной молекулы
 Nmol=ceiling(randn*float(N))
 if (mov==int(Neqv/2.0)) then
    zgran2=KonST*2.0
    zgran1=-KonST
    zdel=zgran2-zgran1
    !print *,'zdel', zdel
    !pause
 endif
 !выбираем что будем делать с молекулой
 steptip=getrand()
 !!запоминаем положение молекулы до
 xdo=x(Nmol) !для центра молекулы
 ydo=y(Nmol)
 zdo=z(Nmol)
 do i=1,Na(Nmol),1
  xado(i)=xa(Nmol,i) !для атома
  yado(i)=ya(Nmol,i)
  zado(i)=za(Nmol,i)
 enddo
 do i=1,bind_kol(Ntip(Nmol)),1              !запоминаем валентные углы
  bind_u_do(i)=bind_u(Nmol,i)
 enddo
 do i=1,koltors(Nmol),1                      !запоминаем торсионные углы
  CosTorsU_do(i)=CosTorsU(Nmol,i)
  TorsU_do(i)=TorsU(Nmol,i)
 enddo

call potenc(Nmol,1)
potenc_do=LJpotenc
Vir_Do=Vir
BindEn_do=Bind_En
TorsEn_do=TorsEnergy
DDLJ_do=DDvir
Tr_en_do=Triple_en
!Перемещение
if (steptip>inter_tip) then
 call tors_ch(Nmol)
 call bind_ch(Nmol)
else
!если рехчастичное
 if (calc_tr==1) then
  if (steptip>0.1) then
   call modif_move(Nmol)
   if (Na(Nmol)>1) then
    call rotate(Nmol)
   endif
  else
   call move(Nmol)
  endif
 else!
  call move(Nmol) !движение
  if (Na(Nmol)>1) then
   call rotate(Nmol)
  endif
 endif
endif
!изменение валентных углов
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call proV(Nmol) !перемещение обратно в объем
call pratom(Nmol)
!нахождение потенциала после
call potenc(Nmol,1)
potenc_posle = LJpotenc
Vir_Posle=Vir
BindEn_posle=Bind_en
TorsEn_posle=TorsEnergy
DDLJ_posle=DDVir
Tr_en_posle=Triple_en
! Нахождение разности потенциалов
deltaP=potenc_posle-potenc_do
deltaVir=Vir_Posle-Vir_Do
deltaBind=BindEn_posle-BindEn_do
deltaTors=TorsEn_posle-TorsEn_do
deltaDDLJ=DDLJ_posle-DDLJ_do
deltaTr=Tr_en_posle-Tr_en_do
!Проверка принятия перемещения
!   if (deltaP>0) then
if (exp(-(deltaP)/Temp)<getrand()) then
 notassept=notassept+1          !переставляем все назад как было
 naDl=naDl+1
 as=0
 !отдельные записи под перемещения
 if (steptip>inter_tip) then
  inter_tot=inter_tot+1
  inter_nac=inter_nac+1
 else
  intra_tot=intra_tot+1
  intra_nac=intra_nac+1
 endif
 !
 x(Nmol)=xdo                            !координаты молекул
 y(Nmol)=ydo
 z(Nmol)=zdo
 do i=1,Na(Nmol),1
  xa(Nmol,i)=xado(i)
  ya(Nmol,i)=yado(i)
  za(Nmol,i)=zado(i)
 enddo
 do i=1,bind_kol(Ntip(Nmol)),1         !валентный угон назад
  bind_u(Nmol,i)=bind_u_do(i)
 enddo
 do i=1,koltors(Nmol),1
  CosTorsU(Nmol,i)=CosTorsU_do(i)
  TorsU(Nmol,i)=TorsU_do(i)
 enddo
 call pratom(Nmol)
else
 !отдельные записи под перемещения
 if (steptip>inter_tip) then
  inter_tot=inter_tot+1
  inter_ac=inter_ac+1
 else
  intra_tot=intra_tot+1
  intra_ac=intra_ac+1
 endif
 !
 assept=assept+1
 aDl=aDl+1
 as=1
 toten=toten+deltap
 TotalVir=TotalVir+deltaVir
 TotBindEn=TotBindEn+deltaBind
 TotTorsEn=TotTorsEn+deltaTors
 TotalDDLJ=TotalDDLJ+deltaDDLJ
 TotalTREn=TotalTrEn+deltaTr
 !call xyzanim()
 call checknear(mov)
endif

if (mov<=Neqv) then !если меньше до равновесия то записываем отдельно
! Neqvil=Neqvil+1
 SampleEqv=SampleEqv+1
 sumP=sumP+toten
 sumP2=sumP2+toten*toten
 sumVir=sumVir+TotalVir
 sumBind=sumBind+TotBindEn
 sumTors=sumTors+TotTorsEn
else if (mov>Neqv) then !записываем на вывод
! Nout=Nout+1
 Sample=Sample+1
 sumPout=sumPout+toten
 sumP2out=sumP2out+(toten+etail1)*(toten+etail1)
 sumVirout=sumVirout+TotalVir
 sumVir2out=sumVir2out+TotalVir*TotalVir
 sumBindout=sumBindout+TotBindEn
 sumTorsout=sumTorsout+TotTorsEn
 sumDDLJ=sumDDLJ+TotalDDLJ
 sumTrEn=sumTrEn+TotalTrEn
endif

if ((mov>Neqv).and.(mod(mov,N*25)==0)) then
 MRDF=MRDF+1
 MRDFN=MRDFN+1
 Mrdf_m=Mrdf_m+1
 call densdistr()
 call bind_ditr()
 !call tors_dist()
 if (calc_rdf==1) then
  call rdf
 endif
 !call rdfm
 if (calc_hend==1) then
  call calc_y
 endif
 call tors_dist()
 !call xyzanim()
endif

if (mod(mov,50000)==0) then
    call xyzanim()
endif

if (mod(mov,Nnew)==0) then
 call outrez
! call calc_y
endif
!!Стартовое изменение шага
!if ((mov<3000).and.(mod(mov,500)==0)) then
! call chDlShag()
!endif
!!Базовое изменениее шага
if (mod(mov,NChan) == 0) then
 call chDlShag()
  !write(6,'(a,i9,a,f9.6,a,f9.6)') 'movie # ', mov, '  max move length: ', DlShag, &
  !&' max rotate angle: ', asin(koefP)*180.0/ChPi
  !write(6,'(4f20.5)') float(mov_ac),float(mov_nac), float(mov_tot),float(mov_ac)/float(mov_tot)
  !write(6,'(4f20.5)') float(rot_ac),float(rot_nac), float(rot_tot),float(rot_ac)/float(rot_tot)
  !write(6,'(4f20.5)') float(int_ac),float(int_nac), float(int_tot),float(int_ac)/float(int_tot)
  if (mov<Neqv) then
   write(6,'(a,f4.2,a,f4.2,a,i7,a,f4.2,a,f4.2,a,i7)') 'intra (ac/nac/tot): ', float(intra_ac)/&
   &float(intra_tot), '/',float(intra_nac)/float(intra_tot),'/', intra_tot, '|  inter (ac/nac&
   &/tot)', float(inter_ac)/float(inter_tot),'/',float(inter_nac)/float(inter_tot),'/', inter_tot
   write(6,'(a,f10.6,a,a,f10.6,a)') 'Max displasement: ', DlShag, trim(razm_rast), '       |  Max &
   &change of bind angle:', asin(BindDeg)*180.0/ChPI, ' º'
   write(6,'(a,f10.6,a,f10.6,a)') 'Max rotate angle: ', asin(RotDeg)*180/ChPI, ' º       |  Max&
   & change of tors angle:', asin(TorsDeg)*180.0/ChPI, ' º'
   write(6,'(a)') '-------------------------------------┴-------------------------------'
  endif
  !!
 inter_tot=0
 inter_ac=0
 inter_nac=0
 intra_tot=0
 intra_nac=0
 intra_ac=0
endif

enddo !конец шага монте карло

!open(13, file='OUTA.txt') !вывод абсолютного положения атомов
!do i=1,N,1
! do j=1,Na(N),1
!  write (13,'(3f30.20)') xaa(i,j),yaa(i,j),zaa(i,j)
! enddo
!enddo
!close(13)
!
!open(13, file='OUTM.txt')
!do i=1,N,1
! write (13,'(3f30.20)') x(i),y(i),z(i)
!enddo
!close(13)
!
!open(13, file='OUTLATM.txt')
! do i=1,N,1
!write (13,'(3f30.20)') x(i)/konst,y(i)/konst,z(i)/konst
!enddo
!close(13)
!!
 do i=1,ckol,1
  do j=i,ckol,1
   write(block_out,'(a,a,a,a)') 'Block-',trim(adjustl(labckol(i))),'-',trim(adjustl&
   &(labckol(j)))
   open(53,file=block_out)
   do k=1,HistRaz,1
    write(53,'(f30.15,$)') Grx(k)
    do bl_no=1,8,1
     write(53,'(a,f30.15,$)') ';',RDF_block(bl_no,(i-1)*ckol+j,k)
    enddo
    write(53,'(a)') ';'
   enddo
   close(53)
   !print *, kc(i,j)
  enddo
 enddo
!!
open(13, file='OUTLATA.xyz') !вывод относительного положения атомов
write(13,'(i3)') totalNa
write(13,'(a52)') 'XYZ file generated by MCP : coordinates in Angstrom'
do i=1,N,1
 do j=1,Na(i),1
  write (13,'(a2,3f20.10)') labela(i,j),(x(i)+xa(i,j))*bs,(y(i)+ya(i,j))*bs,(z(i)+za(i,j))*bs
 enddo
enddo
close(13)
!o2entot=0.0
o2davl=0.0
o2enp=0.0
!нахождение ошибки
do i=1,Mout,1
 o2enp=o2enp+(oenp(Mout)-outsenout)*(oenp(Mout)-outsenout)
 o2davl=o2davl+(odavl(Mout)-outdavlout)*(odavl(Mout)-outdavlout)
enddo
o2enp=sqrt(o2enp/float(Mout))
o2davl=sqrt(o2davl/float(Mout))
!******************************************************************
call cpu_time(finish)
print *, finish-start
!вывод оут файла
call outfile(2)
!запись окончательных знаений в файл
!open(13, file='outy.txt') !переделать вывод
!do i=1,HistRaz,1
! write(13,'(f15.10,a1,f15.10)') Grx(i),';', log(VirGrA(1,i))+(4.0*epsa(1,1,1,1)*(siga(1,1,1,1)**12/Grx(i)**12&
! &-siga(1,1,1,1)**6/Grx(i)**6)/Temp)
!enddo
!close(13)
write(Nazv3,'(a,a1,a)') trim(adjustl(workdir)),'/',trim(adjustl(npr))
 razme=1.0 !be*8.3145107
 razmd=1.0 !1.0/(bs*bs*bs)*be*100.0*1.38065812
open(13, file=Nazv3, position='append')
 write (13,'(7(f25.15,a1),f25.15)') Temp,';', RoNVML,';',outsenout,';',o2enp,';',outdavlout,';',o2davl&
 &,';',outcvout,';',o2cv
close(13)
close(59)
close(56)
!print *, 'chemical potential', -Temp*log(sum_mu(1)/float(nom_mu(1)))
print *, assept
print *, notassept
print *, 'stop'
stop
end program mcprog
!****************************************************************************
subroutine move(Nmol) !передвижение на случайное растояние
use dannie
real(8) randx !случайное передвижение по x
real(8) randy !случайное передвижение по y
real(8) randz !случайное передвижение по z
if ((ptip==3).or.(ptip==5)) then
 if (getrand()>0.5) then
  DlShag=1.0
 else
  DlShag=0.15
 endif
endif
randx=getrand()
randy=getrand()
randz=getrand()
x(Nmol)=x(Nmol)+(randx-0.5)*2*DlShag
y(Nmol)=y(Nmol)+(randy-0.5)*2*DlShag
z(Nmol)=z(Nmol)+(randz-0.5)*2*DlShag
end subroutine move
!------------------------------------------------------------------------------
!процедура создания новой расстановки
subroutine newlatice
use dannie
integer(4) i,j,k !для цикла
integer(4) mi !для молекулы
real(8),allocatable:: sumXm(:) !для нахождения центра масс
real(8),allocatable:: sumYm(:) !для нахождения центра масс
real(8),allocatable:: sumZm(:) !для нахождения центра масс
real(8),allocatable:: sumRm(:) !для нахождения центра масс
real(8),allocatable:: cxm(:) !координата x центра масс
real(8),allocatable:: cym(:) !координата y центра масс
real(8),allocatable:: czm(:) !координата z центра масс
real(8) rastm !максимальное растояние
integer(4),allocatable:: Nmax1(:)
integer(4),allocatable:: Nmax2(:)!номера атомов с максмальным растоянием
real(8) cosfi !cos угла поаорота
real(8) sinfi !sin угла поворота
real(8) cos45 !cos 45
real(8) sin45 !sin 45
real(8) cos54 !cos 54
real(8) sin54 !sin 54
real(8) rast !растояние атома от центра(растояние между двум молекулами)
real(8) xv !временное значение x
real(8) yv !временное значение y
real(8) zv !временное значение z
real(8) xvrem(100,100) !временное значение x (Запомненное)
real(8) yvrem(100,100) !временное значение y
real(8) zvrem(100,100) !временное значение z
real(8) ch_vrx(100,100)
real(8) ch_vry(100,100)
real(8) ch_vrz(100,100)
integer(4) vrem
real(8) maxsproc
integer(4) maxk !индекс вещества с максимальной концентрацией
integer(4) cross !пересечение/непересечение
integer(4) Nshar !шарная молекула
!****************РАССТАНОВКА ЦЕНТРОВ МОЛЕКУЛЫ
allocate(tri_x(NmMax,N*NmMax))     !массив растояний
allocate(tri_y(NmMax,N*NmMax))
allocate(tri_z(NmMax,N*NmMax))
allocate(trx_l(N))
allocate(try_l(N))
allocate(trz_l(N))
allocate(near_mol(N)) !ближайшая молекула
allocate(near_rast(N))
!
select case (tipResh)
case (2) !Гранецентрированная растановка
vrem=1
V=float(N)/RoNV !определяем объем кубика
DlSt=((V**(float(1)/float(3)))/float(storona))
!выделение
allocate(x(N)) !памяти
allocate(y(N)) !под
allocate(z(N)) !массив
!Начало расстановки
do i=0,storona-1,1
 do j=0,storona-1,1
  do k=0,storona-1,1
   !--------------
   x(vrem)=i !временная переменная равна единице все ровно
   y(vrem)=j
   z(vrem)=k
   !--------------
   x(vrem+(storona)**3)=i+0.5
   y(vrem+(storona)**3)=j
   z(vrem+(storona)**3)=k+0.5
   !--------------
   x(vrem+2*(storona)**3)=i
   y(vrem+2*(storona)**3)=j+0.5
   z(vrem+2*(storona)**3)=k+0.5
   !--------------
   x(vrem+3*(storona)**3)=i+0.5
   y(vrem+3*(storona)**3)=j+0.5
   z(vrem+3*(storona)**3)=k
   vrem=vrem+1
  enddo
 enddo
enddo
x=x*DlSt !умножается весь массив
y=y*DlSt
z=z*DlSt
KonSt=V**(float(1)/float(3))
!-------------------------------------------
case (1) !Объемноцентрированная растановка
vrem=1
V=float(N)/RoNV
DlSt=((V**(1.0/3.0))/float(storona))
 !выделение
allocate(x(N)) !памяти
allocate(y(N)) !под
allocate(z(N)) !массив
!Начало расстановки
do i=0,storona-1,1
 do j=0,storona-1,1
  do k=0,storona-1,1
   !---------------
   x(vrem)=i
   y(vrem)=j
   z(vrem)=k
   !---------------
   x(vrem+(storona)**3)=i+0.5
   y(vrem+(storona)**3)=j+0.5
   z(vrem+(storona)**3)=k+0.5
   vrem=vrem+1
  enddo
 enddo
enddo
KonSt=V**(float(1)/float(3))
x=x*DlSt
y=y*DlSt
z=z*DlSt
endselect
!****************************************************************************
!нахождение максимального размера
NmMax=1
NmCh=0    !нахождение максимального
do i=1,Nv,1
 if (Nm(i)>NmMax) then
  NmMax=Nm(i)
 endif
 if (Ntors(i)>NtMax) then
  NtMax=Ntors(i)
 endif
 if (charg_kol(i)>NmCh) then
  NmCh=charg_kol(i)
 endif
enddo
!нахождение центра масс молекулы
allocate(Na(N))
allocate(sumXm(NmMax+1))
allocate(sumYm(NmMax+1))
allocate(sumZm(NmMax+1))
allocate(sumRm(NmMax+1))
allocate(cxm(Nv))
allocate(cym(Nv))
allocate(czm(Nv))
allocate(Nmax1(Nv))
allocate(Nmax2(Nv))
allocate(nom_mu(Nv)) !номер для вычисления хим. потенциала
allocate(sum_mu(Nv)) !сумма для вычисления хим. потенциала
allocate(mu_out(Nv))
allocate(ch_ax(N,NmCh))            !аллокейтим координаты конкретные координаты зарядов
allocate(ch_ay(N,NmCh))            !данные относительно центра
allocate(ch_az(N,NmCh))
allocate(ch_aq(N,NmCh))
allocate(ch_aax(N,NmCh))             !аллоккейтим бсолютные координаты зарядов
allocate(ch_aay(N,NmCh))
allocate(ch_aaz(N,NmCh))
do i=1,Nv,1                      !i-тая молекула j-тый атом
 sumXm(i)=0.0
 sumYm(i)=0.0
 sumZm(i)=0.0
 sumRm(i)=0.0
 Nom_mu(i)=0
 sum_mu(i)=0.0
 do j=1,Nm(i),1
  sumXm(i)=sumXm(i)+xm(i,j)*rm(i,j)
  sumYm(i)=sumYm(i)+ym(i,j)*rm(i,j)
  sumZm(i)=sumZm(i)+zm(i,j)*rm(i,j)
  sumRm(i)=sumRm(i)+rm(i,j)
 enddo
 cxm(i)=sumXm(i)/sumRm(i)
 cym(i)=sumYm(i)/sumRm(i)
 czm(i)=sumZm(i)/sumRm(i)
enddo
!****************************************************************************
!перенос координат атомов в молекуле
do i=1,Nv,1 !i-тая молекула j-тый атом
 do j=1,Nm(i),1               !перенос координат
  xm(i,j)=xm(i,j)-cxm(i)
  ym(i,j)=ym(i,j)-cym(i)
  zm(i,j)=zm(i,j)-czm(i)
 enddo
 do j=1,charg_kol(i),1          !перенос зарядов
  ch_x(i,j)=ch_x(i,j)-cxm(i)
  ch_y(i,j)=ch_y(i,j)-cym(i)
  ch_z(i,j)=ch_z(i,j)-czm(i)
 enddo
enddo

!****************************************************************************
!определение массива координат атомов
!Na=N*NmMax
allocate(Ntip(N))
allocate(koltors(N))
allocate(xa(N,NmMax+1))              !относительные координаты
allocate(ya(N,NmMax+1))
allocate(za(N,NmMax+1))              !абсолютные координаты
allocate(xaa(N,NmMax+1))             !
allocate(yaa(N,NmMax+1))             !
allocate(zaa(N,NmMax+1))             !
allocate(y_rast(N,NmMax+1))          !для расчета y-ка
allocate(CosTorsU(N,NtMax))          !allocate(CosTorsU(N,NtMax)) было раньше
allocate(TorsU(N,NtMax))
!****************************************************************************
!определение самой длинной связи
rast=0.0
rastm=0.0
do i=1,Nv,1
 Nmax2(i)=1 !начальные значения максимального и второго растояния
 Nmax1(i)=1
enddo
!
do mi=1,Nv,1 !mi- тип молекулы(вещество), i-атом
 do i=1,Nm(mi),1
  rast=sqrt(xm(mi,i)**2+ym(mi,i)**2+zm(mi,i)**2)
  if (rast .NE. rastm) then
   Nmax2(mi)=Nmax1(mi)
   Nmax1(mi)=i
  end if
 enddo
!поворот молекулы относительно начала координат так чтобы дальняя молекула была в x 0 0
 if ((xm(mi,Nmax1(mi))/=0).AND.((ym(mi,Nmax1(mi)))/=0)) then
  cosfi=xm(mi,Nmax1(mi))/sqrt(xm(mi,Nmax1(mi))**2+ym(mi,Nmax1(mi))**2) !определение cos угла поворота по z
  sinfi=ym(mi,Nmax1(mi))/sqrt(xm(mi,Nmax1(mi))**2+ym(mi,Nmax1(mi))**2) !определение cos угла поворота по z
!поворот по z
  do i=1,Nm(mi),1
   xv=xm(mi,i)
   yv=ym(mi,i)
   xm(mi,i)= xv*cosfi+yv*sinfi
   ym(mi,i)= (-xv)*sinfi+yv*cosfi
  enddo
  do i=1,Charg_kol(i),1         !Поворот зарядов
   xv=ch_x(mi,i)
   yv=ch_y(mi,i)
   ch_x(mi,i)=xv*cosfi+yv*sinfi
   ch_y(mi,i)=(-xv)*sinfi+yv*cosfi
  enddo
 endif
!
 if ((xm(mi,Nmax1(mi))/=0.0).and.(zm(mi,Nmax1(mi))/=0.0)) then
  cosfi= xm(mi,Nmax1(mi))/sqrt(xm(mi,Nmax1(mi))**2+zm(mi,Nmax1(mi))**2) !определение cos угла поворота по y
  sinfi= zm(mi,Nmax1(mi))/sqrt(xm(mi,Nmax1(mi))**2+zm(mi,Nmax1(mi))**2) !определение sin угла поворота по y
  do i=1,Nm(mi),1
   xv=xm(mi,i)
   zv=zm(mi,i)
   xm(mi,i)= xv*cosfi+zv*sinfi
   zm(mi,i)= (-xv)*sinfi+zv*cosfi
  enddo
  do i=1,Charg_kol(i),1
   xv=ch_x(mi,i)
   zv=ch_z(mi,i)
   ch_x(mi,i)=xv*cosfi+zv*sinfi
   ch_z(mi,i)=(-xv)*sinfi+zv*cosfi
  enddo
 endif
!
 if ((ym(mi,Nmax2(mi))/=0.0) .and. (zm(mi,Nmax2(mi))/=0.0)) then
  cosfi= zm(mi,Nmax2(mi))/sqrt(ym(mi,Nmax2(mi))**2+zm(mi,Nmax2(mi))**2) !определение cos угла поворота по x
  sinfi= ym(mi,Nmax2(mi))/sqrt(ym(mi,Nmax2(mi))**2+zm(mi,Nmax2(mi))**2) !определение sin угла поворота по x
  do i=1,Nm(mi),1
   yv=ym(mi,i)
   zv=zm(mi,i)
   zm(mi,i)= zv*cosfi+yv*sinfi
   ym(mi,i)= (-zv)*sinfi+yv*cosfi
  enddo
  do i=1,charg_kol(i),1
   yv=ch_y(mi,i)
   zv=ch_z(mi,i)
   ch_z(mi,i)=zv*cosfi+yv*sinfi
   ch_y(mi,i)=(-zv)*sinfi+yv*cosfi
  enddo
 endif
enddo !конец для молекулы
!****************************************************************************
!определяем чего больше в смеси
maxsproc=0.0
maxk=1
do i=1,Nv,1
 if (maxsproc<sproc(i)) then
  maxsproc=sproc(i)
  maxk=i
 endif
enddo
!
cos45=sqrt(2.0)/2.0 !косинус 45
sin45=cos45 !синус 45
cos54=1.0/sqrt(3.0) !косинус 54,73
sin54=sqrt(2.0/3.0) !синус 54,73
!запоминаем положение молекул 000
do j=1,Nv,1
 do i=1,Nm(j),1
  xvrem(j,i)=xm(j,i) !запоминаем координаты атомов
  yvrem(j,i)=ym(j,i)
  zvrem(j,i)=zm(j,i)
 enddo
 do i=1,charg_kol(i)            !запоминаем координаты зарядов
  ch_vrx(j,i)=ch_x(j,i)
  ch_vry(j,i)=ch_y(j,i)
  ch_vrz(j,i)=ch_z(j,i)
 enddo
enddo
select case (tipResh)
case (2) !гранецентрированная
!------------------------------------------------------------------
!поворот молекулы по направлению (1,1,1)a
do i=1,Nm(maxk),1 !i-тый атом
 xm(maxk,i)=xvrem(maxk,i)
 ym(maxk,i)=yvrem(maxk,i)
 zm(maxk,i)=zvrem(maxk,i)
 xv=xm(maxk,i)
 yv=ym(maxk,i)
 xm(maxk,i)= xv*cos54-yv*sin54
 ym(maxk,i)= (xv)*sin54+yv*cos54
enddo
do i=1,charg_kol(maxk),1                !вставляем координаты
 ch_x(maxk,i)=ch_vrx(maxk,i)
 ch_y(maxk,i)=ch_vry(maxk,i)
 ch_z(maxk,i)=ch_vrz(maxk,i)
 xv=ch_x(maxk,i)                    !и меняем в соотвествии с углом
 yv=ch_y(maxk,i)
 ch_x(maxk,i)=xv*cos54-yv*sin54
 ch_y(maxk,i)=(xv)*sin54+yv*cos54
enddo
do i=1,Nm(maxk),1
 yv=ym(maxk,i)
 zv=zm(maxk,i)
 zm(maxk,i)= zv*cos45+yv*sin45
 ym(maxk,i)= (-zv)*sin45+yv*cos45
enddo
do i=charg_kol(maxk),1
 yv=ch_vry(maxk,i)
 zv=ch_vrz(maxk,i)
 ch_z(maxk,i)=zv*cos45+yv*sin45
 ch_y(maxk,i)=(-zv)*sin45+yv*cos45
enddo
!распределение атомов (1,1,1) по объему
do i=1,(storona)**3,1
 do j=1,Nm(maxk),1
  xa(i,j)=xm(maxk,j)
  ya(i,j)=ym(maxk,j)
  za(i,j)=zm(maxk,j)
 enddo
 do j=1,charg_kol(maxk),1       !распределение зарядов по 1-му кубу
  ch_ax(i,j)=ch_x(maxk,j)
  ch_ay(i,j)=ch_y(maxk,j)
  ch_az(i,j)=ch_z(maxk,j)
 enddo
enddo
!----------------------------------------------
!распределение атомов (-1,-1,1)d
do i=1,Nm(maxk),1
 xm(maxk,i)=xvrem(maxk,i)
 ym(maxk,i)=yvrem(maxk,i)
 zm(maxk,i)=zvrem(maxk,i)
 xv=xm(maxk,i)
 yv=ym(maxk,i)
 xm(maxk,i)= xv*(-cos54)+yv*(sin54)
 ym(maxk,i)= (-xv)*(sin54)+yv*(-cos54)
enddo
do i=1,charg_kol(maxk),1    !вставляем запомненые координаты
 ch_x(maxk,i)=ch_vrx(maxk,i)
 ch_y(maxk,i)=ch_vry(maxk,i)
 ch_z(maxk,i)=ch_vry(maxk,i)
 xv=ch_x(maxk,i)
 yv=ch_y(maxk,i)
 ch_x(maxk,i)=xv*(-cos54)+yv*(sin54)
 ch_y(maxk,i)=(-xv)*(sin54)+yv*(-cos54)
enddo
do i=1,Nm(maxk),1
 yv=ym(maxk,i)
 zv=zm(maxk,i)
 zm(maxk,i)= zv*cos45-yv*sin45
 ym(maxk,i)= (zv)*sin45+yv*cos45
enddo
do i=1,charg_kol(maxk),1
 yv=ch_y(maxk,i)
 zv=ch_z(maxk,j)
 ch_z(maxk,i)=zv*cos45-yv*sin45
 ch_y(maxk,i)=(zv)*sin45+yv*cos45
enddo
!распределение атомов (1,-1,-1) по объему
do i=1+(storona)**3,2*(storona)**3,1
 do j=1,Nm(maxk),1
  xa(i,j)=xm(maxk,j)
  ya(i,j)=ym(maxk,j)
  za(i,j)=zm(maxk,j)
 enddo
 do j=1,charg_kol(maxk),1
  ch_ax(i,j)=ch_x(maxk,j)
  ch_ay(i,j)=ch_y(maxk,j)
  ch_az(i,j)=ch_z(maxk,j)
 enddo
enddo
!распределение атомов (-1,1,-1)c
do i=1,Nm(maxk),1
 xm(maxk,i)=xvrem(maxk,i)
 ym(maxk,i)=yvrem(maxk,i)
 zm(maxk,i)=zvrem(maxk,i)
 xv=xm(maxk,i)
 yv=ym(maxk,i)
 xm(maxk,i)= xv*(-cos54)-yv*sin54
 ym(maxk,i)= xv*sin54+yv*(-cos54)
enddo
do i=1,charg_kol(maxk),1
 ch_vrx(maxk,i)=ch_x(maxk,i)
 ch_vry(maxk,i)=ch_y(maxk,i)
 ch_vrz(maxk,i)=ch_z(maxk,i)
 xv=ch_x(maxk,i)
 yv=ch_y(maxk,i)
 ch_x(maxk,i)=xv*(-cos54)-yv*sin54
 ch_y(maxk,i)=xv*sin54+yv*(-cos54)
enddo
do i=1,Nm(maxk),1
 yv=ym(maxk,i)
 zv=zm(maxk,i)
 zm(maxk,i)= zv*cos45-yv*sin45
 ym(maxk,i)= (zv)*sin45+yv*cos45
enddo
do i=charg_kol(maxk),1
 yv=ch_y(maxk,i)
 zv=ch_z(maxk,i)
 ch_z(maxk,i)=zv*cos45-yv*sin45
 ch_y(maxk,i)=(zv)*sin45+yv*cos45
enddo
!распределение атомов (-1,1,-1) по объему
do i=1+2*(storona)**3,3*(storona)**3,1
 do j=1,Nm(maxk),1
  xa(i,j)=xm(maxk,j)
  ya(i,j)=ym(maxk,j)
  za(i,j)=zm(maxk,j)
 enddo
 do j=1,charg_kol(maxk),1
  ch_ax(i,j)=ch_x(maxk,j)
  ch_ay(i,j)=ch_y(maxk,j)
  ch_az(i,j)=ch_z(maxk,j)
 enddo
enddo
!распределение атомов (1,-1,-1)b
do i=1,Nm(maxk),1
 xm(maxk,i)=xvrem(maxk,i)
 ym(maxk,i)=yvrem(maxk,i)
 zm(maxk,i)=zvrem(maxk,i)
 xv=xm(maxk,i)
 yv=ym(maxk,i)
 xm(maxk,i)= xv*cos54+yv*sin54
 ym(maxk,i)= -xv*sin54+yv*cos54
enddo
do i=1,charg_kol(maxk),1
 ch_x(maxk,i)=ch_vrx(maxk,i)
 ch_y(maxk,i)=ch_vry(maxk,i)
 ch_z(maxk,i)=ch_vrz(maxk,i)
 xv=ch_x(maxk,i)
 yv=ch_y(maxk,i)
 ch_x(maxk,i)=xv*cos54+yv*sin54
 ch_y(maxk,i)=-xv*sin54+yv*cos54
enddo
do i=1,Nm(maxk),1
 yv=ym(maxk,i)
 zv=zm(maxk,i)
 zm(maxk,i)= zv*cos45+yv*sin45
 ym(maxk,i)= (-zv)*sin45+yv*cos45
enddo
do i=1,charg_kol(maxk),1
 yv=ch_y(maxk,i)
 zv=ch_z(maxk,i)
 ch_z(maxk,i)=zv*cos45+yv*sin45
 ch_y(maxk,i)=(-zv)*sin45+yv*cos45
enddo
!распределение атомов (1,-1,-1) по объему
do i=1+3*(storona)**3,4*(storona)**3,1
 do j=1,Nm(maxk),1
  xa(i,j)=xm(maxk,j)
  ya(i,j)=ym(maxk,j)
  za(i,j)=zm(maxk,j)
 enddo
 do j=1,charg_kol(maxk),1
  ch_ax(i,j)=ch_x(maxk,j)
  ch_ay(i,j)=ch_y(maxk,j)
  ch_az(i,j)=ch_z(maxk,j)
 enddo
enddo
!****************************************************************************
!Объемноцентрированная растановка
case (1)
do i=1,Nm(maxk),1
 xv=xm(maxk,i)
 yv=ym(maxk,i)
 xm(maxk,i)= xv*cos45+yv*sin45
 ym(maxk,i)= (-xv)*sin45+yv*cos45
enddo
do i=1,charg_kol(maxk),1
 xv=ch_x(maxk,i)
 yv=ch_y(maxk,i)
 ch_x(maxk,i)=xv*cos45+yv*sin45
 ch_y(maxk,i)=(-xv)*sin45+yv*cos45
enddo
do i=1,2*(storona)**3,1
 do j=1,Nm(maxk),1
  xa(i,j)=xm(maxk,j)
  ya(i,j)=ym(maxk,j)
  za(i,j)=zm(maxk,j)
 enddo
 do j=1,charg_kol(maxk),1
  ch_ax(i,j)=ch_x(maxk,j)
  ch_ay(i,j)=ch_y(maxk,j)
  ch_az(i,j)=ch_z(maxk,j)
 enddo
enddo
end select
!количество атомов в молекуле под номером
do i=1,N,1
 Na(i)=Nm(maxk)
 Ntip(i)=maxk
enddo
!--------------расстановка остальных молекул---------------
obkol=0
do i=1,Nv,1
 if (i/=maxk) then
  kolNv(i)=int(sproc(i)*N) !количество молекул данного вида
  obkol=obkol+kolNv(i)
 endif
enddo
kolNv(maxk)=N-obkol
do i=1,Nv,1
 sproc(i)=float(kolNv(i))/float(N)
enddo
!
do mi=1,Nv,1            ! рассставляем молекулы
 if (mi/=maxk) then     !которые еще не поставлены
  kpost=0               !обнуляем количество
  do while (kpost<kolNv(mi))    !пока количество поставленных меньше чем необходимо
   randn=getrand()              !выбор случайной молекулы
   Nshar=ceiling(randn*float(N))
!   print *, Ntip(Nshar), randn,sproc(1), sproc(2) ,mi
   !
   if (getrand()<sproc(mi)) then
   if (Ntip(Nshar)==maxk) then
    Ntip(Nshar)=mi
    Na(Nshar)=Nm(mi)
    kpost=kpost+1
    select case (tipResh)
     case (2) !гранецентрированная
     if (Nshar<(storona)**3+1) then
      !поворот молекулы по направлению (1,1,1)a
      do i=1,Nm(mi),1 !i-тый атом
       xm(mi,i)=xvrem(mi,i)
       ym(mi,i)=yvrem(mi,i)
       zm(mi,i)=zvrem(mi,i)
       xv=xm(mi,i)
       yv=ym(mi,i)
       xm(mi,i)= xv*cos54-yv*sin54
       ym(mi,i)= (xv)*sin54+yv*cos54
      enddo
      !поворот заряда
      do i=1,Charg_kol(mi),1    !поворот молекул по направлнию 111
       ch_x(mi,i)=ch_vrx(mi,i)
       ch_y(mi,i)=ch_vry(mi,i)
       ch_z(mi,i)=ch_vrz(mi,i)
       xv=ch_x(mi,i)
       yv=ch_y(mi,i)
       ch_x(mi,i)=xv*cos54-yv*sin54
       ch_y(mi,i)=(xv)*sin54+yv*cos54
      enddo
      do i=1,Nm(mi),1
       yv=ym(mi,i)
       zv=zm(mi,i)
       zm(mi,i)= zv*cos45+yv*sin45
       ym(mi,i)= (-zv)*sin45+yv*cos45
      enddo
      do i=1,charg_kol(mi),1
       yv=ch_y(mi,i)
       zv=ch_z(mi,i)
       ch_z(mi,i)=zv*cos45+yv*sin45
       ch_y(mi,i)=(-zv)*sin45+yv*cos45
      enddo
!распределение атомов (1,1,1) по объему
      do j=1,Nm(mi),1
       xa(Nshar,j)=xm(mi,j)
       ya(Nshar,j)=ym(mi,j)
       za(Nshar,j)=zm(mi,j)
      enddo
      do j=1,charg_kol(mi),1
       ch_ax(Nshar,j)=ch_x(mi,j)
       ch_ay(Nshar,j)=ch_y(mi,j)
       ch_az(Nshar,j)=ch_z(mi,j)
      enddo
     endif
!----------------------------------------------
     if ((Nshar>(storona)**3).and.(Nshar<2*(storona)**3+1)) then
!распределение атомов (-1,-1,1)d
      do i=1,Nm(mi),1
       xm(mi,i)=xvrem(mi,i)
       ym(mi,i)=yvrem(mi,i)
       zm(mi,i)=zvrem(mi,i)
       xv=xm(mi,i)
       yv=ym(mi,i)
       xm(mi,i)= xv*(-cos54)+yv*(sin54)
       ym(mi,i)= (-xv)*(sin54)+yv*(-cos54)
      enddo
      do i=1,charg_kol(mi),1
       ch_x(mi,i)=ch_vrx(mi,i)
       ch_y(mi,i)=ch_vry(mi,i)
       ch_z(mi,i)=ch_vrz(mi,i)
       xv=ch_x(mi,i)
       yv=ch_y(mi,i)
       ch_x(mi,i)=xv*(-cos54)+yv*(sin54)
       ch_y(mi,i)=(-xv)*(sin54)+yv*(-cos54)
      enddo
      do i=1,Nm(mi),1
       yv=ym(mi,i)
       zv=zm(mi,i)
       zm(mi,i)= zv*cos45-yv*sin45
       ym(mi,i)= (zv)*sin45+yv*cos45
      enddo
      do i=1,charg_kol(mi),1
       yv=ch_y(mi,i)
       zv=ch_z(mi,i)
       ch_z(mi,i)=zv*cos45-yv*sin45
       ch_y(mi,i)=(zv)*sin45+yv*cos45
      enddo
!распределение атомов (1,-1,-1) по объему
      do j=1,Nm(mi),1
       xa(Nshar,j)=xm(mi,j)
       ya(Nshar,j)=ym(mi,j)
       za(Nshar,j)=zm(mi,j)
      enddo
      do j=1,charg_kol(mi),1
       ch_ax(Nshar,j)=ch_x(mi,j)
       ch_ay(Nshar,j)=ch_y(mi,j)
       ch_az(Nshar,j)=ch_z(mi,j)
      enddo
     endif
!распределение атомов (-1,1,-1)c
     if ((Nshar>2*(storona)**3).and.(Nshar<3*(storona)**3+1)) then
      do i=1,Nm(mi),1
       xm(mi,i)=xvrem(mi,i)
       ym(mi,i)=yvrem(mi,i)
       zm(mi,i)=zvrem(mi,i)
       xv=xm(mi,i)
       yv=ym(mi,i)
       xm(mi,i)= xv*(-cos54)-yv*sin54
       ym(mi,i)= xv*sin54+yv*(-cos54)
      enddo
      do i=1,charg_kol(mi),1
       ch_x(mi,i)=ch_vrx(mi,i)
       ch_y(mi,i)=ch_vry(mi,i)
       ch_z(mi,i)=ch_vrz(mi,i)
       xv=ch_y(mi,i)
       yv=ch_y(mi,i)
       ch_x(mi,i)=xv*(-cos54)-yv*sin54
       ch_y(mi,i)=xv*sin54+yv*(-cos54)
      enddo
      do i=1,Nm(mi),1
       yv=ym(mi,i)
       zv=zm(mi,i)
       zm(mi,i)= zv*cos45-yv*sin45
       ym(mi,i)= (zv)*sin45+yv*cos45
      enddo
      do i=1,charg_kol(mi),1
       yv=ch_y(mi,i)
       zv=ch_z(mi,i)
       ch_z(mi,i)=zv*cos45-yv*sin45
       ch_y(mi,i)=(zv)*sin45+yv*cos45
      enddo
!распределение атомов (-1,1,-1) по объему
      do j=1,Nm(mi),1
       xa(Nshar,j)=xm(mi,j)
       ya(Nshar,j)=ym(mi,j)
       za(Nshar,j)=zm(mi,j)
      enddo
      do j=1,charg_kol(mi),1
       ch_ax(Nshar,j)=ch_x(mi,j)
       ch_ay(Nshar,j)=ch_y(mi,i)
       ch_az(Nshar,j)=ch_z(mi,i)
      enddo
     endif
     !
     if ((Nshar>3*(storona)**3).and.(Nshar<4*(storona)**3+1)) then
!распределение атомов (1,-1,-1)b
      do i=1,Nm(mi),1
       xm(mi,i)=xvrem(mi,i)
       ym(mi,i)=yvrem(mi,i)
       zm(mi,i)=zvrem(mi,i)
       xv=xm(mi,i)
       yv=ym(mi,i)
       xm(mi,i)= xv*cos54+yv*sin54
       ym(mi,i)= -xv*sin54+yv*cos54
      enddo
      do i=1,charg_kol(mi),1
       ch_x(mi,i)=ch_vrx(mi,i)
       ch_y(mi,i)=ch_vry(mi,i)
       ch_z(mi,i)=ch_vrz(mi,i)
       xv=ch_x(mi,i)
       yv=ch_y(mi,i)
       ch_x(mi,i)=xv*cos54+yv*sin54
       ch_y(mi,i)=-xv*sin54+yv*cos54
      enddo
      do i=1,Nm(mi),1
       yv=ym(mi,i)
       zv=zm(mi,i)
       zm(mi,i)= zv*cos45+yv*sin45
       ym(mi,i)= (-zv)*sin45+yv*cos45
      enddo
      do i=1,charg_kol(mi),1
       yv=ch_y(mi,i)
       zv=ch_z(mi,i)
       ch_z(mi,i)=zv*cos45+yv*sin45
       ch_y(mi,i)=(-zv)*sin45+yv*cos45
      enddo
!распределение атомов (1,-1,-1) по объему
      do j=1,Nm(mi),1
       xa(Nshar,j)=xm(mi,j)
       ya(Nshar,j)=ym(mi,j)
       za(Nshar,j)=zm(mi,j)
      enddo
      do j=1,charg_kol(mi),1
       ch_ax(Nshar,j)=ch_x(mi,i)
       ch_ay(Nshar,j)=ch_y(mi,i)
       ch_az(Nshar,j)=ch_z(mi,i)
      enddo
     endif
!****************************************************************************
!Объемноцентрированная растановка
     case (1)
      do i=1,Nm(mi),1
       xm(mi,i)=xvrem(mi,i)
       ym(mi,i)=yvrem(mi,i)
       zm(mi,i)=zvrem(mi,i)
       xv=xm(mi,i)
       yv=ym(mi,i)
       xm(mi,i)= xv*cos45+yv*sin45
       ym(mi,i)= (-xv)*sin45+yv*cos45
      enddo
      do i=1,charg_kol(mi),1
       ch_x(mi,i)=ch_vrx(mi,i)
       ch_y(mi,i)=ch_vry(mi,i)
       ch_z(mi,i)=ch_vrz(mi,i)
       xv=ch_x(mi,i)
       yv=ch_y(mi,i)
       ch_x(mi,i)=xv*cos45+yv*sin45
       ch_y(mi,i)=(-xv)*sin45+yv*cos45
      enddo
      do j=1,Nm(mi),1
       xa(Nshar,j)=xm(mi,j)
       ya(Nshar,j)=ym(mi,j)
       za(Nshar,j)=zm(mi,j)
      enddo
      do j=1,charg_kol(mi),1        !относительные координаты заряда
       ch_ax(Nshar,j)=ch_x(mi,j)
       ch_ay(Nshar,j)=ch_y(mi,j)
       ch_az(Nshar,j)=ch_z(mi,j)
      enddo
     endselect
     !print *, kpost
    endif
   endif
  enddo
 endif
enddo
!****************************************************************************
!проверка возможности наложения молекул
select case (tipResh)
case (2) !Гранецентрированная
!Если сумма растояний от центра до удаленной молекулы больше межмолекулярного растояния
do i=1,Nv,1
 xv=sqrt((xm(i,Nmax1(i))**2)+(ym(i,Nmax1(i))**2)+(zm(i,Nmax1(i))**2))
 yv=sqrt((xm(i,Nmax2(i))**2)+(ym(i,Nmax2(i))**2)+(zm(i,Nmax2(i))**2))
 zv=(rm(i,Nmax1(i))+rm(i,Nmax2(i)))/2
 if (xv+yv+zv>sqrt(2.0)/2.0*DlSt) then
  cross=1
 endif
enddo
DlShag=sqrt(2.0)/2.0*DlSt/storona !начальная длина шага перемещения
case (1) !Объемноцентрированная
do i=1,Nv,1
 xv=sqrt((xm(i,Nmax1(i))**2)+(ym(i,Nmax1(i))**2)+(zm(i,Nmax1(i))**2))
 yv=sqrt((xm(i,Nmax2(i))**2)+(ym(i,Nmax2(i))**2)+(zm(i,Nmax2(i))**2))
 zv=(rm(i,Nmax1(i))+rm(i,Nmax2(i)))/2
 DlShag=sqrt(3.0)/2.0*DlSt/storona !начальная длинна шага перемещения
 if (xv+yv+zv>DlShag) then
  cross=1
 endif
enddo
end select
!записываем координаты
!open(12, file='NewLatMol.txt')
!do i=1,N,1
!   write (12,'(3f30.15)') x(i),y(i),z(i)
!enddo
!close(12)
!open(13, file='NewLatAtom.txt')
!do i=1,Na,1
!   write (13,'(3f30.15)') xa(i),ya(i),za(i)
!enddo
!close(13)
!****************************************************************************
!абсолютные координаты атомов
totalNa=0
do i=1,N,1
 do j=1,Na(i),1
  totalNa=totalNa+1
  xaa(i,j)=x(i)+xa(i,j)
  yaa(i,j)=y(i)+ya(i,j)
  zaa(i,j)=z(i)+za(i,j)
 enddo
 do j=1,charg_kol(Ntip(i)),1
  ch_aax(i,j)=x(i)+ch_ax(i,j)
  ch_aay(i,j)=y(i)+ch_ay(i,j)
  ch_aaz(i,j)=z(i)+ch_az(i,j)
 enddo
enddo
!
!вставляем этот атом обратно атомы молекулы в объем
!do i=1,N,1
! do j=1,Na(i),1
!  if (xaa(i,j)>KonSt) then
!   xaa(i,j)=xaa(i,j)-KonSt
!  endif
!  if (xaa(i,j)<0) then
!   xaa(i,j)=xaa(i,j)+KonSt
!  endif
!!-------------------------------------
!  if (yaa(i,j)>KonSt) then
!   yaa(i,j)=yaa(i,j)-KonSt
!  endif
!  if (yaa(i,j)<0) then
!   yaa(i,j)=yaa(i,j)+KonSt
!  endif
!!-------------------------------------
!  if (zaa(i,j)>KonSt) then
!   zaa(i,j)=zaa(i,j)-KonSt
!  endif
!  if (zaa(i,j)<0) then
!   zaa(i,j)=zaa(i,j)+KonSt
!  endif
! enddo
!enddo
!--------------------------------------
!записываем абсолютные координаты
!Зписываем название атомов
allocate(labela(N,NmMax+1)) !Буква атома
do i=1,N,1
 do j=1,Na(i)
  labela(i,j)=labelm(Ntip(i),j)
 enddo
enddo
!
open(13, file='NewLatAtomAbs.txt')
do i=1,N,1
 write (13,'(3i4,3f20.15)') i, Ntip(i), Na(i), x(i),y(i),z(i)
 do j=1,Na(i),1
  write (13,'(3f30.15)') xaa(i,j),yaa(i,j),zaa(i,j)
 enddo
 do j=1,charg_kol(Ntip(i)),1
  write(13,'(3f30.15)') ch_aax(i,j),ch_aay(i,j),ch_aaz(i,j)
 enddo
enddo
close(13)

!получли абсолютные координаты атомов
!****************************************************************************
!растановка значений косинусов торсионных углов
do i=1,N,1
 koltors(i)=Ntors(Ntip(i))
 do j=1,Ntors(Ntip(i)),1
  call calc_tors(i,j)
  CosTorsU(i,j)=TorsCosFi !хранятся косинусы углов
  TorsU(i,j)=TorsFi
 enddo
enddo
!определение количества сортов центров
kolcent=0 !общее количество центров во всех атомах
do i=1,Nv,1
 do j=1,Nm(i),1
  kolcent=kolcent+1
 enddo
enddo
allocate(cnames(kolcent))
allocate(nnames(kolcent))
nomer=0
do i=1,Nv,1
 do j=1,Nm(i),1
 nomer=nomer+1
  cnames(nomer)=labelm(i,j) !все лейблы переписаны в один массив
 enddo
enddo

do i=1,kolcent,1
 do j=kolcent,1,-1
  if (cnames(i)==cnames(j)) then
  nnames(i)=j !каждому лейблу соотвествует цифра
  endif
 enddo
enddo
allocate (labckol(kolcent)) !название центра в зависимости от номера
ckol=0
do i=1,kolcent,1
 kolln=0
 do j=1,kolcent,1
  if (i==nnames(j)) then
  kolln=1
  endif
 enddo
 if (kolln==1) then
  ckol=ckol+1
  labckol(ckol)=cnames(i) !
 endif
enddo
!определение количества ФРР
kolfrr=ckol*ckol !размер массива
!определение количества различных центров в молекулах
allocate(Ncrdf(Nv,ckol))
!print *, labckol
do i=1,Nv,1
 do j=1,ckol,1
  Ncrdf(i,j)=0
  do k=1,Nm(i),1
   !print * ,labelm(i,k), labckol(ckol)
   if (labelm(i,k)==labckol(j)) then
    Ncrdf(i,j)=Ncrdf(i,j)+1 !сколько центров типа j  в молекуле i
   endif
  enddo
  !print * , Ncrdf(i,j)
 enddo
enddo
!нахождение доли центров
allocate(Ncsproc(kolfrr))
do i=1,ckol,1
 do k=1,ckol,1
  Ncsum1=0.0
  Ncsum2=0.0
  !sumKNv1=0
  !sumKNv1=0
  rb_kc(i,k)=0.8*sigma(1,1)                         !переделать при необходимости
  do j=1,Nv,1 !по первому проверяемому атому
   Ncsum1=Ncsum1+float(Ncrdf(j,i)*KolNv(j))
   !sumKNv1=sumKNv1+KolNv(j)
  enddo
  do j=1,Nv,1 !по второму проверяемому атому
   Ncsum2=Ncsum2+float(Ncrdf(j,k)*KolNv(j))
   !sumKNv2=sumKNv2+KolNv(j)
  enddo
  Ncsproc((i-1)*ckol+k)=Ncsum1*Ncsum2/float(N)
  !print *, Ncsproc((i-1)*ckol+k)
 enddo
enddo
!углы

allocate(bind_u(N,30))
do i=1,N,1
!print *,bind_kol(Ntip(i))
 do j=1,bind_kol(Ntip(i)),1
  xv1=xa(i,bind1k(Ntip(i),j,1))-xa(i,bind_cent(Ntip(i),j))
  yv1=ya(i,bind1k(Ntip(i),j,1))-ya(i,bind_cent(Ntip(i),j))
  zv1=za(i,bind1k(Ntip(i),j,1))-za(i,bind_cent(Ntip(i),j))
  xv2=xa(i,bind2k(Ntip(i),j,1))-xa(i,bind_cent(Ntip(i),j))
  yv2=ya(i,bind2k(Ntip(i),j,1))-ya(i,bind_cent(Ntip(i),j))
  zv2=za(i,bind2k(Ntip(i),j,1))-za(i,bind_cent(Ntip(i),j))
  r1_bind=sqrt(xv1*xv1+yv1*yv1+zv1*zv1)
  r2_bind=sqrt(xv2*xv2+yv2*yv2+zv2*zv2)
  cos_bind=(xv1*xv2+yv1*yv2+zv1*zv2)/(r1_bind*r2_bind)
  bind_u(i,j)=acos(cos_bind)
 enddo
enddo
!инициализация углов
deg=360.0
kol_deg=360*2
deg_sh=deg/float(kol_deg)
!print *,deg,kol_deg
allocate (hist_u(Nv,10,kol_deg))
allocate (hist_bind(Nv,10,kol_deg))
allocate (hist_norm(Nv,10))
!print* ,'ok'
do i=1,Nv,1 !рассматриваем все молекулы
 do j=1,bind_kol(i),1 !рассматриваем все углы поворота
  hist_norm(i,j)=0
  do k=1,kol_deg,1
   hist_u(i,j,k)=(k-0.5)*deg_sh
   hist_bind(i,j,k)=0
  enddo
 enddo
enddo

!return
end subroutine

subroutine outrez()
use dannie
real(8) eout, dout,outsen,outsdavl
!выводим только после достижения равновесия в файл
if (mov>Neqv) then
 enout=(sumPout)/float(Sample)/float(N)!+LJtail
 viout=(sumVirout)/float(Sample)/float(N)
 davlout=RoNV*Temp-viout*RoNV/3.0!+LJPtail
 bind_out_s=sumBind/float(Sample)/float(N)
 tors_out_s=sumTors/float(sample)/float(N)
 tr_out=sumTrEn/float(sample)/Float(N)
 senout=(senout+enout)
 sdavlout=(sdavlout+davlout)
 Mout=Mout+1
 eout=(enout)*beze+etail1+tr_out !*be*8.3145107
 tr_dout=9.0*tr_out*RoNV
 dout=(davlout)*bezp+ptail1+tr_dout!/(bs*bs*bs)*be*100.0*1.38065812+ptail1
 bind_out_out=bind_out_s*beze
 tors_out_out=tors_out_s*beze
 outsenout=(senout/float(Mout))*beze+etail1!*be*8.3145107+etail1 !средняя потенциальная энергия
 outdavlout=(sdavlout/Mout)*bezp+ptail1!/(bs*bs*bs)*be*100.0*1.38065812+ptail1 !среднее давление
 oenp(Mout)=eout
 odavl(Mout)=dout
 ocv(Mout)=cvout
 cvout=((sumP2out/float(Sample))-(((sumPout+etail1)*(SumPout+etail1))/float(Sample)&
 &/float(Sample)))/Temp/Temp/float(N)*bezcv
 sumCvOut=sumCvOut+CvOut
 OutCvOut=sumCvOut/float(Mout)
 !Вычисление изотермической сжимаемости
 betminus1=-1.0/9.0*(sumVir2out/float(Sample)-sumVirOut/float(Sample)*sumVirOut/float(Sample))&
 &/Temp+1.0/3.0*SumVirOut/float(Sample)+1.0/9.0*SumDDLJ/float(Sample)
 betminus1=betminus1/KonSt/KonSt/KonSt    !в размерном виде не знаю как будет выглядеть
 bette=1/betminus1
 !
 !ошибка
 o2davl=0.0
 o2enp=0.0
 o2cv=0.0
!нахождение ошибки
do i=1,Mout,1
 o2enp=o2enp+(oenp(Mout)-outsenout)*(oenp(Mout)-outsenout)
 o2davl=o2davl+(odavl(Mout)-outdavlout)*(odavl(Mout)-outdavlout)
 o2cv=o2cv+((ocv(Mout)-outcvout)*(ocv(mout)-outcvout))
enddo
o2enp=sqrt(o2enp/float(Mout))
o2davl=sqrt(o2davl/float(Mout))
o2cv=sqrt(o2cv/float(Mout))
 !вывод промежуточной информации
 open(17,file='OUTstep.txt',position='append')
 write (17,'(i10,a1,f30.15,a1,f30.15,a1,f30.15,a1,f30.15,a1,f30.15,a1,f30.15)') Mout*Nnew,&
 &';',eout,';',outsenout,';',o2enp,';',dout,';',outdavlout,';'&
 &,o2davl
 close(17)
 !write (6,'(i8,8(a2,f16.6))') Mout*Nnew,' ',eout,' ',dout,'  ',outsenout,' ',outdavlout, ' '&
 !& ,mu_out(1),' ', mu2_out(1), ' ',DlShag,' ', (float(aDl)+1.0)/(float(naDl)+1.0)
 !write (6,'(a)') '---------------------------------------------------------------------------'
 write (6,'(a, i10,a)') ' step: ', Mout*Nnew, '   ---------------------------------------------&
 &----------------------------'
 write (6,'(a,f15.7,a,a,f15.7,a)') '         energy: ', eout, razm_en,'         pressure: '&
 &, dout, razm_davl
 write (6,'(a,f15.7,a,a,f15.7,a)') ' average energy: ', outsenout, razm_en, ' average pressure: '&
 &, outdavlout, razm_davl
 write (6,'(a,f20.10,a,f20.10)') 'isotermal compressibility: ', bette, ' bT-1  ', betminus1
 write (6, '(a,f9.7,a,f9.7,a,f9.7,a,f9.7)') 'intra:',(eout-bind_out_out-tors_out_out)/eout, ' / inter:'&
 &,(bind_out_out+tors_out_out)/eout,'/ bind : ', bind_out_out/eout, &
 &' / tors: ', tors_out_out/eout
 write (6,'(a,f30.15,a)') 'heat capacity: ', cvout , razm_cv
 write (6,'(a,f30.15,a)') 'average heat capacity: ', OutCvOut, razm_cv
 write (6,'(a, f30.15)') 'dissociation ', float(ndismol/cnkol/nm(1))
 !print *,cvout
 !
   write(6,'(a,f4.2,a,f4.2,a,i7,a,f4.2,a,f4.2,a,i7)') 'intra (ac/nac/tot): ', float(intra_ac)/&
   &float(intra_tot), '/',float(intra_nac)/float(intra_tot),'/', intra_tot, '|  inter (ac/nac&
   &/tot)', float(inter_ac)/float(inter_tot),'/',float(inter_nac)/float(inter_tot),'/', inter_tot
   if ((ptip==3).or.(ptip==5)) then
    write(6,'(a)') 'Warning Max displasement 0.85'
   endif
   write(6,'(a,f10.6,a,a,f10.6,a)') 'Max displasement: ', DlShag, trim(razm_rast), '       |  Max &
   &change of bind angle:', asin(BindDeg)*180.0/ChPI, ' º'
   write(6,'(a,f10.6,a,f10.6,a)') 'Max rotate angle: ', asin(RotDeg)*180/ChPI, ' º       |  Max&
   & change of tors angle:', asin(TorsDeg)*180.0/ChPI, ' º'
   !
 do i=1,ckol,1
  do j=i,ckol,1
   write(6,'(a,a,a,a,a,f20.10)') 'Coordination number ',trim(adjustl(labckol(i))),'-',trim(adjustl&
   &(labckol(j))),' :  ', kc(i,j)
   !!!!перенесение фРРв массив
   do k=1,HistRaz,1
    RDF_block(Mout,(i-1)*ckol+j,k)=VirGrAN((i-1)*ckol+j,k)
    !сразу обнуляем
    SumGrAN((i-1)*ckol+j,k)=0.0
    VirGrAN((i-1)*ckol+j,k)=0.0
   enddo
   !print *, kc(i,j)
  enddo
 enddo
 !
 !!
 do i=1,ckol,1
  do j=i,ckol,1
   write(block_out,'(a,a,a,a)') 'Block-',trim(adjustl(labckol(i))),'-',trim(adjustl&
   &(labckol(j)))
   open(53,file=block_out)
   do k=1,HistRaz,1
    write(53,'(f30.15,$)') Grx(k)
    do bl_no=1,8,1
     write(53,'(a,f30.15,$)') ';',RDF_block(bl_no,(i-1)*ckol+j,k)
    enddo
    write(53,'(a)') ';'
   enddo
   close(53)
   !print *, kc(i,j)
  enddo
 enddo
!!
 !
 if (Calc_hend==1) then
  do i=1,Nv
   write(6,'(a,a,a,f15.10,a)') 'Chemical potential of ', trim(adjustl(MMname(i))),':  '&
   &, mu_out(i), razm_mu
  enddo
 endif
 write(6,'(a,f15.10)') 'accept/not accept : ' , (float(aDl)+1.0)/(float(naDl)+1.0)
 write(6,'(a)') '------------------------------------------------------------------------&
 &-----------------------'
 MRDFN=0
 !
 sumPout=0.0 !начинаем считать заново
 sumP2out=0.0
 sumVirout=0.0
 sumVir2out=0.0
 sumBindout=0.0
 sumTorsout=0.0
 sumDDLj=0.0
 sumTrEn=0.0
 Sample=0
 Nout=0
 !
 call cpu_time(time1)
!
call outfile(1)
!
call cpu_time(time2)

 do i=1,ckol,1 !полная корелляционная функция
  do j=i,ckol,1
   write (rdffile,'(5a)') 'ln(y)-',trim(adjustl(labckol(i))),'-',trim(adjustl(labckol(j))),'.txt'
   write (rdffile2,'(5a)') 'ln(y)-direct-',trim(adjustl(labckol(i))),'-',trim(adjustl(labckol(j))),'.txt'
   open(23, file=rdffile)
   open(24,file=rdffile2)
   do k=1,HistRaz,1 !записываем в функцию
    call potenc_func(i,j,1,1,Grx(k)*Grx(k))
    if (VirGrA((i-1)*ckol+j,k)>0) then
     cavity_out(i,k)=log(VirGrA((i-1)*ckol+j,k))+LJatom/Temp
    else
     cavity_out(i,k)=0
   endif
  enddo
  do k=1,HistRaz,1 !записываем в функцию
   write(24,'(f15.10,a1,f15.10)') Grx(k),';', cavity_out(i,k)
  enddo
  do k=1,ceiling(1.01*maxsig/delta_rast),1 !переносим меньше сигмы
   cavity_out(i,k)=cavity_c_f(i,k)
  enddo
   !cсглаживание выводимой функции
   !3-раза линейное по 3-м точкам
   !5 раз линейное по 5-ти точкам
  do s_k=1,3,1
   smooth(1)=(5.0*cavity_out(i,1)+2.0*cavity_out(i,2)-cavity_out(i,3))/6.0
   smooth(HistRaz)=(5.0*cavity_out(i,HistRaz)+2.0*cavity_out(i,HistRaz-1)&
   &-cavity_out(i,HistRaz-2))/6.0
   do s_i=2,HistRaz-1,1
    smooth(s_i)=(cavity_out(i,s_i-1)+cavity_out(i,s_i)+cavity_out(i,s_i+1))/3.0
   enddo
   do s_i=1,HistRaz,1
    cavity_out(i,s_i)=smooth(s_i)
   enddo
  enddo
   !
  do s_k=1,10,1
   smooth(1)=(3.0*cavity_out(i,1)+2.0*cavity_out(i,2)+cavity_out(i,3)-&
   &cavity_out(i,4))/5.0
   smooth(2)=(4.0*cavity_out(i,1)+3.0*cavity_out(i,2)+2.0*cavity_out(i,3)+&
   &cavity_out(i,4))/10.0
   smooth(HistRaz-1)=(4.0*cavity_out(i,HistRaz)+3.0*cavity_out(i,HistRaz-1)&
   &+2.0*cavity_out(i,HistRaz-2)+cavity_out(i,HistRaz-3))/10.0
   smooth(HistRaz)=(3.0*cavity_out(i,HistRaz)+2.0*cavity_out(i,HistRaz-1)&
   &+cavity_out(i,HistRaz-2)-cavity_out(i,HistRaz-3))/5.0
   do s_i=3,HistRaz-2,1
    smooth(s_i)=(cavity_out(i,s_i-2)+cavity_out(i,s_i-1)+cavity_out(i,s_i)+&
    &cavity_out(i,s_i+1)+cavity_out(i,s_i+2))/5.0
   enddo
   do s_i=ceiling(0.3*maxsig/deltar),HistRaz,1
    cavity_out(i,s_i)=smooth(s_i)
   enddo
  enddo
   !
  do k=1,HistRaz,1 !вывод в файл
    !call potenc_func(i,j,1,1,Grx(k)*Grx(k))
   write(23,'(f15.10,a1,f15.10)') Grx(k),';', cavity_out(i,k)
  enddo
  close(23)
  close(24)
 enddo
enddo
 !
 call calc_cav
 call calc_bridge
 !!вывод информации о торсионном угле
 open(76,file='Tangles.txt')
 write(76,'(a,$)') '#angle'
 do j=1,Nv,1
  do k=1,Ntors(j),1
   write(76,'(a,i1,a,i1,$)') 'Nv',j,'/',k
  enddo
 enddo
 write(76,'(a)') ';'
 do i=1,720,1
  write (76,'(f10.5,a,$)') (i-1)*0.1+0.05,';'
  do j=1,Nv,1
   do k=1,Ntors(j),1
    write(76,'(f20.10,$)') TorsDist(j,k,i)/kolNv(j)/float(MRDF)
   enddo
  enddo
  write(76,'(a)') ';'
 enddo
 close(76)
 !!
 !call saveme
else if (mov<Neqv) then
 en=(sumP)/float(SampleEqv)/float(N)!+LJtail !здесь добавляем хвост!!!!!
 vi=(sumVir)/float(SampleEqv)/float(N)
 davl=RoNV*Temp-vi*RoNV/3.0!+LJPtail
 bind_en_s=(sumBind)/float(SampleEqv)/float(N)
 tors_en_s=sumTors/float(SampleEqv)/float(N)
 sen=(sen+en)
 sdavl=(sdavl+davl)
 sbind=sbind+bind_en_s
 stors=stors+tors_en_s
 Meqvil=Meqvil+1
 eout=(en)*beze+etail1 !*be*8.3145107
 dout=(davl)*bezp+ptail1!/(bs*bs*bs)*be*100.0*1.38065812
 bind_en_out=bind_en_s*beze
 tors_en_out=tors_en_s*beze
 outsen=(sen/Meqvil)*beze+etail1!*be*8.3145107
 outsdavl=(sdavl/Meqvil)*bezp+ptail1!/(bs*bs*bs)*be*100.0*1.38065812
 cv=((sumP2/float(SampleEqv))-(sumP*sumP)/float(SampleEqv)/float(SampleEqv))/Temp/Temp/N*bezcv
 !write (6,'(i10,4(a2,f20.10))') Meqvil*Nnew,' ',eout,' ',dout, ' ',outsen,' ',outsdavl
 write (6,'(a, i10,a)') ' step: ', Meqvil*Nnew, '   ---------------------------------------------'
 write (6,'(a,f15.7,a,a,f15.7,a)') '         energy: ', eout, razm_en,'         pressure: '&
 &, dout, razm_davl
 write (6,'(a,f15.7,a,a,f15.7,a)') ' average energy: ', outsen, razm_en, ' average pressure: '&
 &, outdavl, razm_davl
 write (6,'(a,f20.10,a)') 'heat capacity: ', cv , razm_cv
 write (6, '(f15.8,a,f15.8,a,f15.8)')eout, '/', bind_en_out, &
 &' ', tors_en_out
 print *, eout, bind_en_out, tors_en_out
 !
   write(6,'(a,f4.2,a,f4.2,a,i7,a,f4.2,a,f4.2,a,i7)') 'intra (ac/nac/tot): ', float(intra_ac)/&
   &float(intra_tot), '/',float(intra_nac)/float(intra_tot),'/', intra_tot, '|  inter (ac/nac&
   &/tot)', float(inter_ac)/float(inter_tot),'/',float(inter_nac)/float(inter_tot),'/', inter_tot
   if ((ptip==3).or.(ptip==5)) then
    write(6,'(a)') 'Warning Max displasement 0.85'
   endif
   write(6,'(a,f10.6,a,a,f10.6,a)') 'Max displasement: ', DlShag, trim(razm_rast), '       |  Max &
   &change of bind angle:', asin(BindDeg)*180.0/ChPI, ' º'
   write(6,'(a,f10.6,a,f10.6,a)') 'Max rotate angle: ', asin(RotDeg)*180/ChPI, ' º       |  Max&
   & change of tors angle:', asin(TorsDeg)*180.0/ChPI, ' º'
 !
 write(6,'(a,f15.10)') 'accept/not accept : ' , (float(aDl)+1.0)/(float(naDl)+1.0)
 write(6,'(a)') '-------------------------------------------------------------------------&
 &--------------------------'
 !print *,outsen,'<<',outsdavl
 !write (6,'(i10,f30.20)') Meqvil*Nnew,sumP!, sumVir, float(Nnew)
 !начинаем заново считать
 sumP=0.0
 SumP2=0.0
 sumPout=0.0
 sumVir=0.0
 sumTors=0.0
 sumBind=0.0
 Neqvil=0
 SampleEqv=0
else if (mov==Neqv) then
 write (6,'(a70)') '****************************start productation*********************&
 &************************************'
endif
!call saveme
end subroutine

subroutine rdf()
use dannie
integer(4) i,j,ii,jj
integer(4) bin,bin2
real(8) vHist !число гистограммы
real(8) x1,x2,y1,y2,z1,z2
real(8) rastHista
!вывод на экран
!*******************************Молекулярная
!обнуляем гистограмму
do i=1,HistRaz,1
 Hist(i)=0.0
 Gr(i)=0.0
enddo
do i=1,22,1
 Hist_w(i)=0.0
 Gr_w(i)=0.0
enddo
!находим гистограмму
do i=1,N,1
 do j=i,N,1
  if (i/=j) then
   call rastmol(i,j)
   !определяем только на половине длинны объема
   if (rastmm<rGrcut*rGrcut) then
    rastmm=sqrt(rastmm)
    bin=ceiling((rastmm-0.5*deltar)/deltar)
    Hist(bin)=Hist(bin)+2.0
    !проверка для молекул с ямами
    !if (rastmm>rxk1_w) then
    ! print *,rastmm,rxk1_w
    ! pause
    ! if (rastmm<rxk2_w) then
    !  print *,rastmm,rxk2_w
    !  pause
    ! endif
    !endif
    if ((rastmm>rxk1_w).and.(rastmm<rxk2_w)) then       !rxk1_w,rxk2_w
     bin2=ceiling((rastmm-rxk1_w)/dr_well)
     !print *, bin, rastmm,(rastmm-rxk1_w)/dr_well, dr_well, rxk2_w
     !pause
     Hist_w(bin2)=Hist_w(bin2)+2.0
    endif
   endif
  endif
 enddo
enddo
!нахождение функции распределения из гистограммы
do i=1,HistRaz,1
 Gr(i)=Hist(i)/Nideal(i)/float(N)
 SumGr(i)=SumGr(i)+Gr(i)
 VirGr(i)=SumGr(i)/float(Mrdf)
enddo
    !Эта фишка используется только для модели центральных сил
do i=1,122,1                                    !поэтому находится только здесь
 Gr_w(i)=hist_w(i)/Nideal_well(i)/float(N)
 SumGr_w(i)=SumGr_w(i)+Gr_w(i)
 VirGr_w(i)=SumGr_w(i)/float(Mrdf)
enddo
open(12,file='rdf_w.txt')
do i=1,122,1
 write(12,'(f15.10,a1,f15.10)') Grx_well(i),';',VirGr_w(i)
enddo
close(12)
open(12, file='RDFM.txt')
do i=1,HistRaz,1
 write(12,'(f15.10,a1,f15.10)') Grx(i),';',VirGr(i)
enddo
close(12)
!**************Атомная
do inom=1,ckol,1
 do jnom=inom,ckol,1
  !обнуляем гистограмму
  do i=1,HistRaz,1
   HistAtom((inom-1)*ckol+jnom,i)=0.0
   GrA((inom-1)*ckol+jnom,i)=0.0
  enddo
  !находим гистограмму
  do i=1,N,1
   do ii=1,Na(i),1 !номер проверяемого центра в молекуле i
    do j=i,N,1 !номер молекулы с которой проверяем
     if (i/=j) then !проверяем только центры в разных молекулах
      do jj=1,Na(j),1 !номер центра в молекуле с которой проверяем
       !квадрат растояния между двумя атомами
       if ((labelm(Ntip(i),ii)==labckol(inom)).and.(labelm(Ntip(j),jj)==labckol(jnom))) then
        x1=xaa(i,ii)
        x2=xaa(j,jj)
        y1=yaa(i,ii)
        y2=yaa(j,jj)
        z1=zaa(i,ii)
        z2=zaa(j,jj)
        !___________________    минимальное растояние
        if (x2-x1>KonSt/2) then
         x2=x2-KonSt
        else if (x1-x2>KonSt/2) then
         x2=x2+KonSt
        endif
        !___________________
        if (y2-y1>KonSt/2) then
         y2=y2-KonSt
        else if (y1-y2>KonSt/2) then
         y2=y2+KonSt
        endif
        !___________________
        if (z2-z1>KonSt/2) then
         z2=z2-KonSt
        else if (z1-z2>KonSt/2) then
         z2=z2+KonSt
        endif
        !___________________
        dax=(x2-x1)
        day=(y2-y1)
        daz=(z2-z1)
        !___________________
        rastHista=(dax*dax+day*day+daz*daz) !квадрат растояния
        if (rastHista<rGrcut*rGrcut) then
         bin=ceiling((sqrt(rastHista)-0.5*deltar)/deltar)
         HistAtom((inom-1)*ckol+jnom,bin)=HistAtom((inom-1)*ckol+jnom,bin)+2.0
        endif
       endif
      enddo
     endif
    enddo
   enddo
  enddo
 enddo
enddo
!нахождение функции распределения из гистограммы
do inom=1,ckol,1
 do jnom=inom,ckol,1
  do i=1,HistRaz,1
   vHist=HistAtom((inom-1)*ckol+jnom,i)
   GrA((inom-1)*ckol+jnom,i)=vHist/(Nideal(i))/Ncsproc((inom-1)*ckol+jnom)
   SumGrA((inom-1)*ckol+jnom,i)=SumGrA((inom-1)*ckol+jnom,i)+GrA((inom-1)*ckol+jnom,i)
   VirGrA((inom-1)*ckol+jnom,i)=SumGrA((inom-1)*ckol+jnom,i)/float(Mrdf)
   SumGrAN((inom-1)*ckol+jnom,i)=SumGrAN((inom-1)*ckol+jnom,i)+GrA((inom-1)*ckol+jnom,i)
   VirGrAN((inom-1)*ckol+jnom,i)=SumGrAN((inom-1)*ckol+jnom,i)/float(MrdfN)
  enddo
 enddo
enddo
!вывод средней гистограмы на печать
!сглаживание ФРР
do inom=1,ckol,1
 do jnom=inom,ckol,1
  kc(inom,jnom)=0.0     !обнуляем координационное число
  write (rdffile,'(a4,a,a1,a,a4)') 'RDF-',trim(adjustl(labckol(inom))),'-',trim(adjustl(labckol(jnom))),'.txt'
  open(23, file=rdffile)
  do i=1,HistRaz,1
   write(23,'(f15.10,a1,f15.10)') Grx(i),';', VirGrA((inom-1)*ckol+jnom,i)
   if (Grx(i)<Rb_kc(inom,jnom)) then
    kc(inom,jnom)=kc(inom,jnom)+(VirGrA((inom-1)*ckol+jnom,i+1)+VirGrA((inom-1)*ckol+jnom,i))&
    &/(Grx(i+1)-Grx(i))/2.0*(Grx(i+1)+Grx(i))*(Grx(i+1)+Grx(i))
   endif
  enddo
  kc(inom,jnom)=kc(inom,jnom)*ChPI*RoNV*Ncsproc((inom-1)*ckol+jnom)
  close(23)
 enddo
enddo
do i=1,ckol,1 !полная корелляционная функция
 do j=i,ckol,1
  write (rdffile,'(5a)') 'H-',trim(adjustl(labckol(i))),'-',trim(adjustl(labckol(j))),'.txt'
  !write (rdffile2,'(5a)') 'H-2-',trim(adjustl(labckol(i))),'-',trim(adjustl(labckol(j))),'.txt'
  open(23, file=rdffile)
  !open(24,file=rdffile2)
  do k=1,HistRaz,1
   total_c_f((i-1)*ckol+j,k)=VirGrA((i-1)*ckol+j,k)-1.0
   write(23,'(f15.10,a1,f15.10)') Grx(k),';', total_c_f((i-1)*ckol+j,k)
   !if (mod(k,2)==0) then
   ! write(24,'(f15.10,a1,f15.10)') Grx(k),';', total_c_f((i-1)*ckol+j,k)
   !endif
  enddo
  close(23)
  !close(24)
 enddo
enddo
call smooth_rdf()
end subroutine
!!Модификация нахождения функции распреления
!
subroutine potenc(Nmol,fromN) !нахождение потенциала от одной молекулы
use dannie !используем уже размеченные переменные
integer(4) pi,pj,fromN
integer(4) LJmi,LJmj
integer(4) rai

LJpotenc=0.0
Vir=0.0
DDVir=0.0
min_rast_m=KonSt*KonSt*100.0
!call pratom(Nmol) !получили все атомы в объеме !необходимо в дальнейшем для РДФ
do pj=1,Na(Nmol),1
 xprov1(pj)=xa(Nmol,pj)+x(Nmol)
 yprov1(pj)=ya(Nmol,pj)+y(Nmol)
 zprov1(pj)=za(Nmol,pj)+z(Nmol)
enddo !переносим данные проверяемой молекулы
!
do pi=fromN,N,1
 if (pi .NE. Nmol) then
  !call rastmol(Nmol,pi)
  if (x(Nmol)-x(pi) > KonSt/2) then
   xv=x(pi)+KonSt
   trx_l(pi)=KonSt    !этот максимум находится для одной молекулы
  else if (x(pi)-x(Nmol) > KonSt/2) then
   xv=x(pi)-KonSt
   trx_l(pi)=-KonSt
  else
   xv=x(pi)
   trx_l(pi)=0.0
  endif
!print *, 'ololo1'
  if (y(Nmol)-y(pi) > KonSt/2) then
   yv=y(pi)+KonSt
   try_l(pi)=KonSt
  else if (y(pi)-y(Nmol) > KonSt/2) then
   yv=y(pi)-KonSt
   try_l(pi)=-KonSt
  else
   yv=y(pi)
   try_l(pi)=0.0
  endif
!print *, 'ololo2'
  if (z(Nmol)-z(pi) > zdel/2) then
   zv=z(pi)+zdel
   trz_l(pi)=zdel
  else if (z(pi)-z(Nmol) > zdel/2) then
   zv=z(pi)-zdel
   trz_l(pi)=-zdel
  else
   zv=z(pi)
   trz_l(pi)=0.0
  endif
!рассчитываем растояния
!do rai=1,Na(pi),1
! tr_dx=()-()
! Drast(Nmol*NmMax+pj,pi*NmMax+rai)=1 !решить что тут
!enddo
!виртуальные координаты молекулы
  do rai=1,Na(pi),1
   xprov2(rai)=xa(pi,rai)+xv
   yprov2(rai)=ya(pi,rai)+yv
   zprov2(rai)=za(pi,rai)+zv
  enddo
!растояние
  dmx=(-x(Nmol)+xv)
  dmy=(-y(Nmol)+yv)
  dmz=(-z(Nmol)+zv)
  rastmm=(dmx*dmx+dmy*dmy+dmz*dmz)
  if (rastmm<min_rast_m) then
   min_rast_m=rastmm
   near_mol(Nmol)=pi
  endif
  if (rastmm .LT. RadCut2) then
   LJmol=0.0
   VirLJmol=0.0
   DDmol=0.0
   do LJmi=1,Na(Nmol),1
    do LJmj=1,Na(pi),1
    !print *, 'ololo4'
     dax=(-xprov1(LJmi)+xprov2(LJmj))
     day=(-yprov1(LJmi)+yprov2(LJmj))
     daz=(-zprov1(LJmi)+zprov2(LJmj))
     rastaa=(dax*dax+day*day+daz*daz)   !вот же находим растояние
     !надо сюда прикрутить
     !tri_x(LJmi,pi*LJmj)=dax        !растояние записываем в массив
     !tri_y(LJmi,pi*LJmj)=day
     !tri_z(LJmi,pi*LJmj)=daz
!     print *, 'ok'
!    pause
     provKof=1.0/(rastaa)*(dax*dmx+day*dmy+daz*dmz)
     call potenc_func(Ntip(Nmol),Ntip(pi),LJmi,LJmj,rastaa)
     LJmol=LJmol+LJatom
     VirLJmol=VirLJmol+VirLJatom*ProvKof
     DDmol=DDmol+DDLJatom
    enddo
   enddo
   LJpotenc=LJpotenc+LJmol !Нахождение потенциала между атомами 2х моелкул
   Vir=Vir+VirLJmol
   DDVir=DDvir+DDmol
  endif
 endif
enddo
DbEn=LJpotenc
DbVir=Vir
call bind_energy(Nmol)
LJpotenc=LJpotenc+bind_en
call tors_energy(Nmol)
LJpotenc=LJpotenc+torsEnergy
!print *, 'ololo5'
if (calc_tr==1) then
 call triple_pot(Nmol)
 Triple_en=Uijk
 LJpotenc=LJpotenc+Triple_en
endif
end subroutine
!
!Подпрограмма проверяет на нахождение центра молекулы в заданном объеме
subroutine proV(Nmol)
use dannie
if (x(Nmol)>KonSt) then
 x(Nmol)=x(Nmol)-KonSt
endif
if (x(Nmol)<0) then
 x(Nmol)=x(Nmol)+KonSt
endif
!-------------------------------------
if (y(Nmol)>KonSt) then
 y(Nmol)=y(Nmol)-KonSt
endif
if (y(Nmol)<0) then
 y(Nmol)=y(Nmol)+KonSt
endif
!-------------------------------------
if (z(Nmol)>zgran2) then
 z(Nmol)=z(Nmol)-zdel
endif
if (z(Nmol)<zgran1) then
 z(Nmol)=z(Nmol)+zdel
endif
!--------------------------------------
end subroutine

!*****************************************для атомной
!начальная растановка атомов
subroutine creatAA()
use dannie
integer(4) i,j
!находим измененные абсолютные коордиинаты атома
!абсолютные координаты атомов
do i=1,N,1 !номер молекулы
 do j=1,Na(i),1 !номер атома в молекуле
!вставляем этот атом обратно атомы молекулы в объем
  if (xaa(i,j)>KonSt) then
   xaa(i,j)=xaa(i,j)-KonSt
  endif
  if (xaa(i,j)<0) then
   xaa(i,j)=xaa(i,j)+KonSt
  endif
!-------------------------------------
  if (yaa(i,j)>KonSt) then
   yaa(i,j)=yaa(i,j)-KonSt
  endif
  if (yaa(i,j)<0) then
   yaa(i,j)=yaa(i,j)+KonSt
  endif
!-------------------------------------
  if (zaa(i,j)>zgran2) then
   zaa(i,j)=zaa(i,j)-zdel
  endif
  if (zaa(i,j)<zgran1) then
   zaa(i,j)=zaa(i,j)+zdel
  endif
 enddo
enddo
end subroutine
!***************************************************************
!проверка атома на нахождение в объеме
subroutine pratom(Nmolp)
use dannie
integer(4) i, Nmolp
do i=1,Na(Nmolp),1
!номер атома для проверки
 xaa(Nmolp,i)=x(Nmolp)+xa(Nmolp,i)
 yaa(Nmolp,i)=y(Nmolp)+ya(Nmolp,i)
 zaa(Nmolp,i)=z(Nmolp)+za(Nmolp,i)
!вставляем этот атом обратно атомы молекулы в объем
 if (xaa(Nmolp,i)>KonSt) then
  xaa(Nmolp,i)=xaa(Nmolp,i)-KonSt
 endif
 if (xaa(Nmolp,i)<0) then
  xaa(Nmolp,i)=xaa(Nmolp,i)+KonSt
 endif
!-------------------------------------
 if (yaa(Nmolp,i)>KonSt) then
  yaa(Nmolp,i)=yaa(Nmolp,i)-KonSt
 endif
 if (yaa(Nmolp,i)<0) then
  yaa(Nmolp,i)=yaa(Nmolp,i)+KonSt
 endif
!-------------------------------------
 if (zaa(Nmolp,i)>zgran2) then
  zaa(Nmolp,i)=zaa(Nmolp,i)-zdel
 endif
 if (zaa(Nmolp,i)<zgran1) then
  zaa(Nmolp,i)=zaa(Nmolp,i)+zdel
 endif
!--------------------------------------
enddo
end subroutine
!Функция рандомирует случайные числа
!начало рандома
subroutine randomn
use generator
integer(4) i
real(8) sseed1,sseed2,sseed3
integer(4) seed(20)
character(20) a,b,c
real(8) n11,n21,n31
call date_and_time(a,b,c,seed)
sseed1=0.0
sseed2=0.0
sseed3=0.0
!print *,seed
do i=1,20,1 !здесь шарной бред
 sseed1=sseed1+seed(i)
 if (mod(i,2)==0) then
  sseed2=sseed2+seed(i)
 else
  sseed2=sseed2/2.0
 endif
 if (mod(i,3)==0) then
  sseed3=sseed3+seed(i)
 else
  sseed3=sseed3-seed(i)/2.2
 endif
enddo
!print *, sseed1,sseed2,sseed3
!pause
!проверка  тригонометрических функций
n11=abs(cos(sseed1/1000))
n21=abs(sin(sseed2/1000))
n31=abs(cos(sseed3/1000))
!call random_seed(int(sseed))
!call random_number(n11)
!call random_number(n21)
!call random_number(n31)
m1=2147483563
a1=40014
b1=12345
m2=2147483399
a2=40692
b2=54321
max3=2147483647
max1=m1-1
max2=m2-1
n1=ceiling((n11*max1)-1)
n2=ceiling((n21*max2)-1)
n3=ceiling((n31*max3)-1)
!   open(15,file='rand.txt')
 do i=1,500
 call randstart()
 randmass(i)=outrand
!       write(15,'(f30.20)') randmass(i)
 enddo
!write (6,'(a10,3i15)') 'zatravka',n1,n2,n3
!   close(15)
end subroutine
!*******************************************************
!Рандомная функция начало
subroutine  randstart()
use generator
n1=abs(mod(a1*n1+b1,m1))
n2=abs(mod(a2*n2+b2,m2))
n3=abs(n3*1664525+1013904223)
if (float(n3)/float(max3)<0.5) then
 outrand=float(n1)/float(max1)
else
 outrand=float(n2)/float(max2)
endif
end subroutine
!*******************************************************
!Рандомная функция
function getrand()
use generator
integer(4) repick
!n1=abs(mod(a1*n1+b1,m1))
!n2=abs(mod(a2*n2+b2,m2))
n3=abs(n3*1664525+1013904223)
repick=ceiling((float(n3)/float(max3))*500)
getrand=randmass(repick)
call randstart()
randmass(repick)=outrand
!getrand=0.012
end function
!нахождение минимадьного растояния между двумя атомами
!
subroutine rastatom(rasti,rastj)
use dannie
integer(4) rasti,rastj
 dax=(-xprov1(rasti)+xprov2(rastj))
 day=(-yprov1(rasti)+yprov2(rastj))
 daz=(-zprov1(rasti)+zprov2(rastj))
 rastaa=(dax*dax+day*day+daz*daz)
end subroutine
!****************************************************************************
!Функция нахождения растояния между молекулами
subroutine rastmol(Nmol1,Nmol2)
use dannie
integer(4) rai
!real(8) rast
if (x(Nmol1)-x(Nmol2) > KonSt/2) then
 xv=x(Nmol2)+KonSt
else if (x(Nmol2)-x(Nmol1) > KonSt/2) then
 xv=x(Nmol2)-KonSt
else
 xv=x(Nmol2)
endif

if (y(Nmol1)-y(Nmol2) > KonSt/2) then
 yv=y(Nmol2)+KonSt
else if (y(Nmol2)-y(Nmol1) > KonSt/2) then
 yv=y(Nmol2)-KonSt
else
 yv=y(Nmol2)
endif

if (z(Nmol1)-z(Nmol2) > zdel/2.0) then
 zv=z(Nmol2)+zdel
else if (z(Nmol2)-z(Nmol1) > zdel/2.0) then
 zv=z(Nmol2)-zdel
else
 zv=z(Nmol2)
endif
!виртуальные координаты молекулы
do rai=1,Na(Nmol2),1
 xprov2(rai)=xa(Nmol2,rai)+xv
 yprov2(rai)=ya(Nmol2,rai)+yv
 zprov2(rai)=za(Nmol2,rai)+zv
enddo
!растояние
dmx=(-x(Nmol1)+xv)
dmy=(-y(Nmol1)+yv)
dmz=(-z(Nmol1)+zv)
rastmm=(dmx*dmx+dmy*dmy+dmz*dmz)
end subroutine
!****************************************************************************
!Поворот на случайный угол
!****************************************************************************
subroutine rotate(Nmol) !поворот на случаный угол
use dannie
real(8) randalfa !случайный поворот относительно оси
real(8) randbeta !случайный поворот относительно оси
real(8) randgamma !случайный поворот относительно оси
real(8) sinOx,sinOy,sinOz
real(8) cosOx,cosOy,cosOz
integer(4) Nmol,i

io=0
randalfa=getrand()
randbeta=getrand()
randgamma=getrand()
!Находим синусы и косинусы поворачиваемых углов
sinOx=RotDeg*(randalfa-0.5)*2.0
cosOx=sqrt(1.0-sinOx*sinOx)
sinOy=RotDeg*(randbeta-0.5)*2.0
cosOy=sqrt(1.0-sinOy*sinOy)
sinOz=RotDeg*(randgamma-0.5)*2.0
cosOz=sqrt(1.0-sinOz*sinOz)
do i=1,Na(Nmol),1
 xv=xa(Nmol,i)
 yv=ya(Nmol,i)
 zv=za(Nmol,i)
 xa(Nmol,i)=cosOx*cosOy*xv+(cosOx*sinOy*sinOz-sinOx*cosOz)*yv+(cosOx*sinOy*cosOz+sinOx*sinOz)*zv
 ya(Nmol,i)=sinOx*cosOy*xv+(sinOx*sinOy*sinOz+cosOx*cosOz)*yv+(sinOx*sinOy*cosOz-cosOx*sinOz)*zv
 za(Nmol,i)=-sinOy*xv+cosOy*sinOz*yv+cosOy*cosOz*zv
enddo

!!Поворачиваем молекулу относительно оси x
!do i=1,Na(Nmol),1
! yv=ya(Nmol,i)
! zv=za(Nmol,i)
! ya(Nmol,i)= yv*cosOx+zv*sinOx
! za(Nmol,i)= (-yv)*sinOx+zv*cosOx
!enddo
!do i=1,charg_kol(Ntip(Nmol)),1     !!поворот зарядов
! yv=ch_ay(Nmol,i)
! zv=ch_az(Nmol,i)
! ch_ay(Nmol,i)=yv*cosOx+zv*sinOx
! ch_az(Nmol,i)=(-yv)*sinOx+zv*cosOx
!enddo
!!Поворачиваем молекулу относительно оси y
!do i=1,Na(Nmol),1
! xv=xa(Nmol,i)
! zv=za(Nmol,i)
! xa(Nmol,i)= xv*cosOy+zv*sinOy
! za(Nmol,i)= (-xv)*sinOy+zv*cosOy
!enddo
!do i=1,charg_kol(Ntip(Nmol)),1     !!поворот заряда
! xv=ch_ax(Nmol,i)
! zv=ch_az(Nmol,i)
! ch_ax(Nmol,i)=xv*cosOy+zv*sinOy
! ch_az(Nmol,i)=(-xv)*sinOy+zv*cosOy
!enddo
!!Поворачиваем молекулу относительно оси z
!do i=1,Na(Nmol),1
! xv=xa(Nmol,i)
! yv=ya(Nmol,i)
! xa(Nmol,i)= xv*cosOz+yv*sinOz
! ya(Nmol,i)= (-xv)*sinOz+yv*cosOz
!enddo
!do i=1,charg_kol(Ntip(Nmol)),1
! xv=ch_ax(Nmol,i)
! yv=ch_ay(Nmol,i)
! ch_ax(Nmol,i)=xv*cosOz+yv*sinOz
! ch_ay(Nmol,i)=(-xv)*sinOz+yv*cosOz
!enddo
end subroutine rotate
!****************************************************************************
subroutine tails()
use dannie
integer(4) i,j,mi,mj
real(8) sig3
real(8),allocatable:: xa1(:,:)
real(8),allocatable:: ya1(:,:)
real(8),allocatable:: za1(:,:)
real(8),allocatable:: xa2(:,:)
real(8),allocatable:: ya2(:,:)
real(8),allocatable:: za2(:,:)
real(8),allocatable:: outrtail(:)
real(8),allocatable:: outnrot(:)
real(8),allocatable:: outrot(:)
real(8),allocatable:: outcos(:)
real(8),allocatable:: outtaile(:)
real(8),allocatable:: outnrote(:)
real(8) rc, Rcut2,Rcut2m
real(8) randalfa !случайный поворот относительно оси
real(8) randbeta !случайный поворот относительно оси
real(8) randgamma !случайный поворот относительно оси
real(8) NRtail,RotTail,sumTail,NRtailE
real(8) sintx,costx,sinty,costy,sintz,costz
real(8)  sumcos, cosg !,ttail2,
integer(4) kstep, step
real(8) dxa,dya,dza,dxm,dym,dzm
cosg=0.0
LJtail=0.0
LJPtail=0.0

allocate (xa1(Nv,NmMax))
allocate (ya1(Nv,NmMax))
allocate (za1(Nv,NmMax))
allocate (xa2(Nv,NmMax))
allocate (ya2(Nv,NmMax))
allocate (za2(Nv,NmMax))
allocate (outrtail(int(20/drtail)))
allocate (outtaile(int(20/drtail)))
!allocate (outnrot(int(20/drtail)))
!allocate (outnrote(int(20/drtail)))
allocate (outrot(int(20/drtail)))
allocate (outcos(int(20/drtail)))
do mi=1,Nv,1 !инициализация начальных координат
 do i=1,Nm(mi),1
  xa1(mi,i)=xm(mi,i)
  ya1(mi,i)=ym(mi,i)
  za1(mi,i)=zm(mi,i)
  xa2(mi,i)=xm(mi,i)
  ya2(mi,i)=ym(mi,i)
  za2(mi,i)=zm(mi,i)
 enddo
enddo
! введение монтекарло для решения интеграла
Rcut2=RadCut
do kstep=1,int(20/drtail),1 !относительная разница усреднения по растоянию
! подсчет без кручения
!NRtail=0.0
!NRtailE=0.0
!do mi=1,Nv,1
! do mj=1,Nv,1
!   do i=1,Nm(mi),1
!    do j=1,Nm(mj),1
!     sig3=siga(mi,mj,i,j)
!     sig3=sig3*sig3*sig3
!     rc=sig3/Rcut2/Rcut2/Rcut2 !третья степень
!     rc=rc*rc !шестая степень
!     NRtailE=Rcut2*Rcut2*epsa(mi,mj,i,j)*(rc*rc-rc)*sproc(mi)*sproc(mj)+NRtailE
!     NRtail=Rcut2*Rcut2*epsa(mi,mj,i,j)*(rc*rc-0.5*rc)*sproc(mi)*sproc(mj)+NRtail
!    enddo
!   enddo
!   !добавить процент
! enddo
!enddo
!NRtailE=NRtailE*2.0*ChPI*RoNV*4.0
!NRtail=NRtail*2.0/3.0*ChPI*48.0*RoNV*RoNV
!подсчет с кручением
!расчет среднего методом монте-карло
sumTail=0.0
sumTailE=0.0
sumcos=0.0
dxm=Rcut2 ! не надо вычислять несколько раз
dym=0.0
dzm=0.0
do mi=1,Nv,1
 do mj=1,Nv,1
  do step=1,Nstep,1
   !поворачиваем одну молекулу
   randalfa=getrand()
   randbeta=getrand()
   randgamma=getrand()
   !Находим синусы и косинусы поворачиваемых углов
   sintx=(randalfa-0.5)*2
   costx=sqrt(1-sintx*sintx)
   sinty=(randbeta-0.5)*2
   costy=sqrt(1-sinty*sinty)
   sintz=(randgamma-0.5)*2
   costz=sqrt(1-sintz*sintz)
   !Поворачиваем молекулу относительно оси x
   do i=1,Nm(mi),1
    yv=ya1(mi,i)
    zv=za1(mi,i)
    ya1(mi,i)= yv*costx+zv*sintx
    za1(mi,i)= (-yv)*sintx+zv*costx
   enddo
   !Поворачиваем молекулу относительно оси y
   do i=1,Nm(mi),1
    xv=xa1(mi,i)
    zv=za1(mi,i)
    xa1(mi,i)= xv*costx+zv*sintx
    za1(mi,i)= (-xv)*sintx+zv*costx
   enddo
   !Поворачиваем молекулу относительно оси z
   do i=1,Nm(mi),1
    xv=xa1(mi,i)
    yv=ya1(mi,i)
    xa1(mi,i)= xv*costx+yv*sintx
    ya1(mi,i)= (-xv)*sintx+yv*costx
   enddo
   !поворачиваем другую молекулу
   randalfa=getrand()
   randbeta=getrand()
   randgamma=getrand()
   !Находим синусы и косинусы поворачиваемых углов
   sintx=(randalfa-0.5)*2
   costx=sqrt(1-sintx*sintx)
   sinty=(randbeta-0.5)*2
   costy=sqrt(1-sinty*sinty)
   sintz=(randgamma-0.5)*2
   costz=sqrt(1-sintz*sintz)
  !Поворачиваем молекулу относительно оси x
   do i=1,Nm(mj),1
    yv=ya2(mj,i)
    zv=za2(mj,i)
    ya2(mj,i)= yv*costx+zv*sintx
    za2(mj,i)= (-yv)*sintx+zv*costx
   enddo
   !Поворачиваем молекулу относительно оси y
   do i=1,Nm(mj),1
    xv=xa2(mj,i)
    zv=za2(mj,i)
    xa2(mj,i)= xv*costx+zv*sintx
    za2(mj,i)= (-xv)*sintx+zv*costx
   enddo
   !Поворачиваем молекулу относительно оси z
   do i=1,Nm(mj),1
    xv=xa2(mj,i)
    yv=ya2(mj,i)
    xa2(mj,i)= xv*costx+yv*sintx
    ya2(mj,i)= (-xv)*sintx+yv*costx
   enddo
   !нахождение с учетом углов
   !
   !print *,xa1(1),ya1(1),za1(1)
    RotTail=0.0
    RotTailEn=0.0
    do i=1,Nm(mi),1
     do j=1,Nm(mj),1
      dxa=xa2(mj,j)-xa1(mi,i)+Rcut2
      dya=ya2(mj,j)-ya1(mi,i)
      dza=za2(mj,j)-za1(mi,i)
      Rcut2m=dxa*dxa+dya*dya+dza*dza
      cosg=(dxa*dxm+dya*dym+dza*dzm)/Rcut2m
      call potenc_func(Ntip(mi),Ntip(mj),i,j,Rcut2m)
      !print *, Rcut2m,mi,mj,i,j, LJatom
      RotTailEn=RotTailEn+Rcut2*Rcut2*Sproc(mi)*Sproc(mj)*LJatom
      RotTail=RotTail+Rcut2*Rcut2*Sproc(mi)*Sproc(mj)*VirLJatom*cosg
     enddo
    enddo
    sumTail=sumTail+RotTail
    sumcos=sumcos+cosg
    sumtailE=sumtailE+RotTailEn
    !
   enddo
   !конец степа
  enddo
enddo
outcos(kstep)=sumcos/float(Nstep)
outrtail(kstep)=Rcut2
outtaile(kstep)=sumtaile/float(Nstep)*2.0*ChPI*RoNV
outrot(kstep)=sumTail/float(Nstep)*2.0/3.0*ChPI*RoNV*RoNV
Rcut2=Rcut2+drtail
enddo
!сравнение
ttail=0.0
ttaile=0.0
do i=1,int(20/drtail),1
 if ((i/=1) .or. (i/=int(20/drtail))) then
  ttail=ttail+outrot(i)
  ttaile=ttaile+outtaile(i)
 else
  ttail=ttail+outrot(i)*0.5
  ttaile=ttaile+outtaile(i)*0.5
 endif
enddo
!нахождение 3-ей части давления !для леннард джонса не работает - нет аналитического решения
if (ptip==2) then
 etail2=0.0
 ptail2=0.0
 do mi=1,Nv,1
  do mj=1,Nv,1
   do i=1, Nm(mi),1
    do j=1,Nm(mj),1
     sig3=siga(mi,mj,i,j) !Было для ЛД
     sig3=sig3*sig3*sig3
     rc=sig3/(RadCut+20.0)/(RadCut+20.0)/(RadCut+20.0)
     etail2=epsa(mi,mj,i,j)*sproc(mi)*sproc(mj)*sig3*(rc*rc*rc/3.0-rc)+etail2
     ptail2=epsa(mi,mj,i,j)*sproc(mi)*sproc(mj)*sig3*(2.0/3.0*rc*rc*rc-rc)+ptail2
    enddo
   enddo
  enddo
 enddo
 etail2=etail2*8.0/3.0*ChPI*RoNV!*be*8.3145107 !размерный вид потом
 ptail2=ptail2*16.0/3.0*ChPI*RoNV*RoNV!/(bs*bs*bs)*be*100.0*1.38065812
 etail1=etail1+etail2
 ptail1=ptail1+ptail2
endif
etail1=ttaile*drtail*beze!*be*8.3145107 !поставить/убрать размерный вид
ptail1=-ttail*drtail*bezp!/(bs*bs*bs)*be*100.0*1.38065812
write (6,'(a,f20.10,a)') 'tail correction for energy ', etail1, razm_en
write (6,'(a,f20.10,a)') 'tail correction for pressure ',ptail1, razm_davl

end subroutine

subroutine totalen()
use dannie
integer(4) i
real(8) rr
TotEn=0.0
TotDbEn=0.0
TotBindEn=0.0
TotTorsEn=0.0
TotalVir=0.0
TotalDDlj=0.0
TotalTrEn=0.0
do i=1,N,1
 call potenc(i,1)
 TotDbEn=TotDbEn+DbEn
 TotalVir=TotalVir+DbVir
 TotBindEn=TotBindEn+Bind_En
 TotTorsEn=TotTorsEn+TorsEnergy
 TotalDDLJ=TotalDDLJ+DDVir
 TotalTrEn=TotalTrEn+Triple_en
enddo
TotDBEn=TotDbEn/2.0
TotalVir=TotalVir/2.0
TotalTrEn=TotalTrEn/3.0
TotEn=TotDBEn+TotalTrEn+TotBindEn+TotTorsEn
!print *, toten, totbinden, tottorsen
!print *,toten
write (6,'(a,f30.10,a)') 'Total energy of start latice per molecule: ', toten*beze/&
&float(N), razm_en
write (6,'(a,f30.10,a)') 'Total triple energy of start latice per molecule: ', TotalTrEn*beze/&
&float(N), razm_en
print *, TotalTrEn
open(19,file='potenc')
do i=1,5000,1
 rr=(float(i)*0.001)*(float(i)*0.001)
 call potenc_func(1,1,1,1,rr)
 if (LJatom<100.0) then
  write(19,'(f20.10,a1,f20.10,a1,f20.10,a1,f20.10,a,f20.10)') sqrt(rr),';',LJatom,';',VirLJatom,&
  & ';', VirLJatom/sqrt(rr), ';', DDLJatom
 endif
enddo
close(19)
prov_var=1
end subroutine

subroutine tors_ch(Nmol) !поворот относительно торсионого угла
use dannie
integer(4) ito,jto,kto
real(8) xv,yv,zv,delx,dely,delz
real(8) cos1,cos2,sin1,sin2
real(8) CosPov,SinPov

do ito=1,Ntors(Ntip(Nmol)),1 !делаем пока не окнчатся торсионные связи
!!проверка
!!call calc_tors(Nmol,ito)
!print *,'do',TorsCosFi,TorsFi
!!
 delx=xa(Nmol,fsta(Ntip(Nmol),ito))
 dely=ya(Nmol,fsta(Ntip(Nmol),ito))
 delz=za(Nmol,fsta(Ntip(Nmol),ito))
 do jto=1,Na(Nmol),1 !переставляем первый атом в начало координат
  xa(Nmol,jto)=xa(Nmol,jto)-delx
  ya(Nmol,jto)=ya(Nmol,jto)-dely
  za(Nmol,jto)=za(Nmol,jto)-delz
 enddo
 !переставляем вторую молекулу на ось
 sin1=ya(Nmol,seca(Ntip(Nmol),ito))/sqrt(ya(Nmol,seca(Ntip(Nmol),ito))*ya(Nmol,seca(Ntip(Nmol),ito))&
 &+za(Nmol,seca(Ntip(Nmol),ito))*za(Nmol,seca(Ntip(Nmol),ito)))
 cos1=za(Nmol,seca(Ntip(Nmol),ito))/sqrt(ya(Nmol,seca(Ntip(Nmol),ito))*ya(Nmol,seca(Ntip(Nmol),ito))&
 &+za(Nmol,seca(Ntip(Nmol),ito))*za(Nmol,seca(Ntip(Nmol),ito)))
 do jto=1,Na(Nmol),1
  yv=-za(Nmol,jto)*sin1+ya(Nmol,jto)*cos1
  zv=za(Nmol,jto)*cos1+ya(Nmol,jto)*sin1
  za(Nmol,jto)=zv
  ya(Nmol,jto)=yv
 enddo
 sin2=za(Nmol,seca(Ntip(Nmol),ito))/sqrt(xa(Nmol,seca(Ntip(Nmol),ito))*xa(Nmol,seca(Ntip(Nmol),ito))&
 &+za(Nmol,seca(Ntip(Nmol),ito))*za(Nmol,seca(Ntip(Nmol),ito)))
 cos2=xa(Nmol,seca(Ntip(Nmol),ito))/sqrt(xa(Nmol,seca(Ntip(Nmol),ito))*xa(Nmol,seca(Ntip(Nmol),ito))&
 &+za(Nmol,seca(Ntip(Nmol),ito))*za(Nmol,seca(Ntip(Nmol),ito)))
 do jto=1,Na(Nmol),1
  xv=xa(Nmol,jto)*cos2+za(Nmol,jto)*sin2
  zv=-xa(Nmol,jto)*sin2+za(Nmol,jto)*cos2
  za(Nmol,jto)=zv
  xa(Nmol,jto)=xv
 enddo
if (getrand()<0.5) then
  SinPov=(getrand()-0.5)*2.0*TorsDeg !синус угла на который нужно поворачивать
   CosPov=sqrt(1.0-SinPov*SinPov) !косинус угла
 else
  SinPov=0.86602540378443900000
  CosPov=-0.5
endif
 CosIsh=CosTorsU(Nmol,ito) !косинус угла который был до поворота
 SinIsh=sqrt(1.0-CosIsh*CosIsh) !синус угла до поворота
 CosTorsU(Nmol,ito)=CosPov*CosIsh-SinPov*SinIsh
 !поворот части молекулы
 do jto=1,nrota(Ntip(N),ito),1
  yv=ya(Nmol,rota(Ntip(Nmol),ito,jto))*CosPov+za(Nmol,rota(Ntip(Nmol),ito,jto))*SinPov
  zv=-za(Nmol,rota(Ntip(Nmol),ito,jto))*CosPov+ya(Nmol,rota(Ntip(Nmol),ito,jto))*SinPov
  ya(Nmol,rota(Ntip(Nmol),ito,jto))=yv
  za(Nmol,rota(Ntip(Nmol),ito,jto))=zv
 enddo
 !обратный поворот по одной осе
 do jto=1,Na(Nmol),1
  xv=xa(Nmol,jto)*cos2-za(Nmol,jto)*sin2
  zv=xa(Nmol,jto)*sin2+za(Nmol,jto)*cos2
  za(Nmol,jto)=zv
  xa(Nmol,jto)=xv
 enddo
 !обратный поворот по другой оси
 do jto=1,Na(Nmol),1
  yv=ya(Nmol,jto)*Cos1+za(Nmol,jto)*sin1
  zv=-ya(Nmol,jto)*sin1+za(Nmol,jto)*cos1
  ya(Nmol,jto)=yv
  za(Nmol,jto)=zv
 enddo
 do jto=1,Na(Nmol),1 !переставляем атомы обратно как было
  xa(Nmol,jto)=xa(Nmol,jto)+delx
  ya(Nmol,jto)=ya(Nmol,jto)+dely
  za(Nmol,jto)=za(Nmol,jto)+delz
 enddo
!!проверка
!call calc_tors(Nmol,ito)
!print *,'posle',TorsCosFi,TorsFi
!print *,'---'
!!
enddo

end subroutine

subroutine tors_energy(Nom)
use dannie
integer tui,tuj,Nom
real(8) coso,coso2,coso3,coso4,coso5
torsEnergy=0.0
do tui=1,koltors(Nom),1
 call calc_tors(Nom,tui)   !вычисляем торсионный угол и его косинус
 CosTorsU(Nom,tui)=TorsCosFi
 TorsU(Nom,tui)=TorsFi
 !print *,CosTorsU(Nom,tui),TorsU(Nom,tui)
 coso=CosTorsU(Nom,tui)
 coso2=2.0*coso*coso-1.0
 coso3=4.0*coso*coso*coso-3.0*coso
 coso4=8.0*coso*coso*coso*coso-8.0*coso*coso+1.0
 coso5=16.0*coso*coso*coso*coso*coso-20.0*coso*coso*coso+5.0*coso
 torsEnergy=torsEnergy+TorsKoef(Ntip(Nom),tui,1)*(1.0+coso)+TorsKoef(Ntip(Nom),tui,2)*(1.0-coso2)+&
 &TorsKoef(Ntip(Nom),tui,3)*(1.0+coso3) !+TorsKoef(Ntip(Nom),tui,4)*coso4+&
 !&TorsKoef(Ntip(Nom),tui,5)*coso5
 !print *,CosTorsU(Nom,tui),TorsU(Nom,tui),torsEnergy
 !pause
enddo
end subroutine

subroutine calc_y()
use dannie
real(8) x1,y1,z1,sum_energy
real(8) vxran,vyran,vzran,rast_vrem,xv,yv,zv
real(8) sin_a,cos_a,sin_b,cos_b
integer(4) bin,i,j,k,mu_n
integer(4) rast
real(8) mu2_potenc
!print *,vxran,vyran,vzran
do k=1,Nv,1
 do mu_n=1,1000,1
  vxran=getrand()*KonSt
  vyran=getrand()*KonSt
  vzran=getrand()*KonSt
  sum_energy=0.0
  do i=1,N,1
   x1=x(i)!+xa(i,j)
   y1=y(i)!+ya(i,j)
   z1=z(i)!+za(i,j)
   !___________________ минимальное растояние
   if (vxran-x1>KonSt/2.0) then
    x1=x1+KonSt
   else if (x1-vxran>KonSt/2.0) then
    x1=x1-KonSt
   endif
   !___________________
   if (vyran-y1>KonSt/2.0) then
    y1=y1+KonSt
   else if (y1-vyran>KonSt/2.0) then
    y1=y1-KonSt
   endif
   !___________________
   if (vzran-z1>KonSt/2.0) then
    z1=z1+KonSt
   else if (z1-vzran>KonSt/2.0) then
    z1=z1-KonSt
   endif
   !___________________
   dax=(vxran-x1)
   day=(vyran-y1)
   daz=(vzran-z1)
   rast_vrem=dax*dax+day*day+daz*daz
   Y_rast(i,1)=rast_vrem !квадрат растояния
   if (rast_vrem < RadCut2) then
    call potenc_func(Ntip(i),k,1,1,rast_vrem)
    sum_energy=sum_energy+LJatom
   endif
  enddo
  !print *, sum_energy
  nom_mu(k)=nom_mu(k)+1
  sum_mu(k)=sum_mu(k)+exp(-(sum_energy+etail1/float(N))/Temp)
  !do i=1,N,1
  ! if (k==Ntip(i)) then
  !  if (Y_rast(i,1)< RadCut2) then
  !   !print *, histraz, bin
  !   bin=ceiling(sqrt(y_rast(i,1))/deltar)
  !   nom_y(k,bin)=nom_y(k,bin)+1
  !   call potenc_func(Ntip(i),k,1,1,y_rast(i,1))
  !   potenc_minus=sum_energy-LJatom
  !   HistY(k,bin)=HistY(k,bin)+exp(-(potenc_minus+etail2/float(N))/Temp)
  !   sum_minus(k,bin)=sum_minus(k,bin)+potenc_minus
  !  endif
  ! endif
  !enddo
 enddo
 mu_out(k)=-Temp*log(sum_mu(k)/float(nom_mu(k)))
  !print *, mu_out(k), sum_mu(k), nom_mu(k), temp
  !write (yfile,'(a,a,a)') 'y-hist-',trim(adjustl(labckol(k))),'.txt'
  !open(29, file=yfile)
  ! do i=1,HistRaz,1
  !  write(29,'(f15.10,a1,f60.30)') Grx(i),';',log(HistY(k,i)/float(nom_y(k,i)))+mu_out(k)/Temp
  ! enddo
  !close(29)
enddo

!другой алгоритм
do k=1,Nv,1
 do mu_n=1,100,1
 Nk_rand=ceiling(getrand()*N)
 do while (Ntip(Nk_rand)/=k) !выбираем случайную молекулу данного типа
  Nk_rand=ceiling(getrand()*N)
 enddo
 y2_nom=y2_nom+1
 do rast=1,ceiling(1.05*maxsig/delta_rast),1 ! !kon_rast,1!
  sin_a=(getrand()-0.5)*2.0!поворот соединяющего вектора
  cos_a=sqrt(1.0-sin_a*sin_a)
  sin_b=(getrand()-0.5)*2.0
  cos_b=sqrt(1.0-sin_b*sin_b)
  vyran=0.0 !координаты по y и z остаются такими же
  vzran=0.0
  vxran=delta_rast*float(rast)
  xv=vxran !поворот
  zv=vzran
  vxran= xv*cos_a+zv*sin_a
  vzran=(-xv)*sin_a+zv*cos_a
  xv=vxran
  yv=vyran
  vxran= xv*cos_b+yv*sin_b
  vyran=(-xv)*sin_b+yv*cos_b
  vxran=vxran+x(Nk_rand)
  vyran=vyran+y(Nk_rand)
  vzran=vzran+z(Nk_rand)
  y_potenc2(k,rast)=0.0 !обнуляем потенциал
  mu2_potenc=0.0
  if (vxran>KonSt) then
   vxran=vxran-KonSt
  endif
  if (vyran>KonSt) then
   vyran=vyran-KonSt
  endif
  if (vzran>KonSt) then
   vzran=vzran-KonSt
  endif!_____
  if (vxran<0.0) then
   vxran=vxran+KonSt
  endif
  if (vyran<0.0) then
   vyran=vyran+KonSt
  endif
  if (vzran<0.0) then
   vzran=vzran+KonSt
  endif
  do i=1,N,1
   !if (i/=Nk_rand) then !если молекула не та же самая
    x1=x(i)!+xa(i,j)
    y1=y(i)!+ya(i,j)
    z1=z(i)!+za(i,j)
    !___________________    минимальное растояние
    if (vxran-x1>KonSt/2.0) then
     x1=x1+KonSt
    else if (x1-vxran>KonSt/2.0) then
     x1=x1-KonSt
    endif
    !___________________
    if (vyran-y1>KonSt/2.0) then
     y1=y1+KonSt
    else if (y1-vyran>KonSt/2.0) then
     y1=y1-KonSt
    endif
    !___________________
    if (vzran-z1>KonSt/2.0) then
     z1=z1+KonSt
    else if (z1-vzran>KonSt/2.0) then
     z1=z1-KonSt
    endif
    !___________________
    dax=(vxran-x1)
    day=(vyran-y1)
    daz=(vzran-z1)
    rast_vrem=dax*dax+day*day+daz*daz
    if (rast_vrem<RadCut2) then !если входит в радиус то считем потенциал
     call potenc_func(Ntip(i),k,1,1,rast_vrem)
     if (i/=Nk_rand) then
      y_potenc2(k,rast)=y_potenc2(k,rast)+LJatom
     endif
     !mu2_potenc=mu2_potenc+LJatom
    endif
   !endif
  enddo
  !mu2_nom(k)=mu2_nom(k)+1
  !mu2_sum(k)=mu2_sum(k)+exp(-(mu2_potenc+etail1/float(N))/Temp)
  sum_potenc2(k,rast)=sum_potenc2(k,rast)+exp((-y_potenc2(k,rast)+etail1/float(N))/Temp)
 enddo
 enddo !mu
 !mu2_out(k)=-Temp*log(mu2_sum(k)/float(mu2_nom(k)))
 !print *, mu2_out(k),mu_out(k),mu2_nom(k),mu2_sum(k)
 !вывод информации
 do i=1,ceiling(1.05*maxsig/delta_rast),1
  cavity_c_f(k,i)=log(sum_potenc2(k,i)/float(y2_nom))+mu_out(k)/Temp
 enddo
  write (yfile,'(a,a,a)') 'y2-',trim(adjustl(labckol(k))),'.txt'
 open(29, file=yfile)
  do i=1,ceiling(1.05*maxsig/delta_rast),1!!kon_rast,1 !  !
   write(29,'(f20.10,a1,f30.15)') Grx(i),';', cavity_c_f(k,i)
  enddo
 close(29)
 !smooth
 do sm_i=1,10,1
  smooth(1)=(5.0*cavity_c_f(k,1)+2.0*cavity_c_f(k,2)-cavity_c_f(k,3))/6.0
  smooth(ceiling(1.05*maxsig/delta_rast))=(5.0*cavity_c_f(k,ceiling(1.05*&
  &maxsig/delta_rast))+2.0*cavity_c_f(k,ceiling(1.05*maxsig/delta_rast)-1)&
  &-cavity_c_f(k,ceiling(1.5*maxsig/delta_rast)-2))/6.0
  do i=2,ceiling(1.05*maxsig/delta_rast)-1,1
   smooth(i)=(cavity_c_f(k,i-1)+cavity_c_f(k,i)+cavity_c_f(k,i+1))/3.0
  enddo
  !do i=1,ceiling(1.05*maxsig/delta_rast),1
  ! cavity_c_f(k,i)=smooth(i)
  !enddo
 enddo

 !write (yfile,'(a,a,a)') 'y2-2-',trim(adjustl(labckol(k))),'.txt'
 !open(29, file=yfile)
 ! do i=1,kon_rast,1
 !  if (Mod(i,2)==1) then
 !   !cavity_c_f(k,i)=log(sum_potenc2(k,i)/float(y2_nom))+mu2_out(k)/Temp
 !   write(29,'(f20.10,a1,f30.15)') i*delta_rast,';', cavity_c_f(k,i)
 !  endif
 ! enddo
 !close(29)
enddo
end subroutine

subroutine potenc_func(N1,N2,i1,i2,RastSQR)
use dannie
integer(4) N1,N2,i1,i2
real(8) RastSQR
real(8) rp,wp,kp,lr,ep,sp
real(8) du1,ddu1,dul,ddul,du2,ddu2
LJatom=0.0 !зануляем потенциал
VirLJatom=0.0
if (ptip==1) then !Карра-Коновалова
 rz=sqrt(RastSQR)
 rz2=1.0/RastSQR
 rz6=rz2*rz2*rz2
 rz7=rz6/rz
 expp=exp(alf(N1,N2,i1,i2)*(1.0-rz/siga(N1,N2,i1,i2)))
 LJatom=(kk1(N1,N2,i1,i2)*expp+kk2(N1,N2,i1,i2))*rz6
 VirLJatom=(dkk1(N1,N2,i1,i2)*expp*rz6+dkk2(N1,N2,i1,i2)*expp*rz7+&
 &dkk3(N1,N2,i1,i2)*rz7+dkk4(N1,N2,i1,i2)*rz7)*rz
 !print *,LJatom,RastSQR,rz,expp,rz6,rz7
 !pause
 !---------------
endif
if (ptip==2) then !Леннарда Джонса
 rz2=siga(N1,N2,i1,i2)*siga(N1,N2,i1,i2)/RastSQR
 rz6=rz2*rz2*rz2
 LJatom=4.0*epsa(N1,N2,i1,i2)*(rz6*rz6-rz6)
 VirLJatom=epsa(N1,N2,i1,i2)*(-rz6*rz6 + 0.5*rz6)*48.0
endif
if (ptip==3) then !леннард джонс с ямой
 rz=sqrt(RastSQR)
 rp=0.75 !r_pot(Ntip(N1),Ntip(N2),i1,i2)
 wp=0.015 !w_pot(Ntip(N1),Ntip(N2),i1,i2)
 kp=20000.0 !k_pot(Ntip(N1),Ntip(N2),i1,i2)
 !lp=lrazm !l_pot(Ntip(N1),Ntip(N2),i1,i2)
 lr=0.5*(1.0+tanh((rz-rp)/wp))
 ep=epsa(N1,N2,i1,i2)
 sp=siga(N1,N2,i1,i2)
 rz2=sp*sp/RastSQR
 rz6=rz2*rz2*rz2
 Uinter=4.0*ep*(rz6*rz6-rz6)
 Uintra=0.5*kp*(rz-lp)*(rz-lp)+Dep
 LJatom=(1.0-lr)*Uintra+lr*Uinter
 sep=0.5/wp/((cosh((rz-rp)/wp))*(cosh((rz-rp)/wp)))
 VirLJatom=kp*(rz-lp)*(1-lr)+4.0*ep*(6.0*rz6-12.0*rz6*rz6)/rz*lr&
 &+Uintra*sep+4.0*ep*(rz6*rz6-rz6)*sep
 VirLJatom=VirLjatom*rz
 du2=4.0*ep/rz*(6.0*rz6-12.0*rz6*rz6)
 ddu2=4.0*ep/rz/rz*(156.0*rz6*rz6-42.0*rz6)
 dul=0.5/wp/((cosh((rz-rp)/wp))*(cosh((rz-rp)/wp)))
 ddul=-tanh((rz-rp)/wp)/wp/wp/((cosh((rz-rp)/wp))*(cosh((rz-rp)/wp)))
 du1=kp*(rz-lp)
 ddu1=kp
 !!
 VirLJatom=lr*du2+(1-lr)*du1+Uinter*dul-Uintra*dul
 VirLJatom=VirLJatom*rz
 !выражение для производной
 DDLJatom=rz*rz*lr*ddu2+(2.0*rz*rz*dul+rz*lr)*du2+(rz*rz-rz*rz*lr)*ddu1+&
 &(-2.0*rz*rz*dul-rz*lr+rz)*du1+(rz*rz*Uinter-rz*rz*Uintra)*ddul+&
 &(rz*Uinter-rz*Uintra)*dul
 !print *,LJatom,RastSQR,rz,expp,rz6,rz7
endif
if (ptip==4) then
rz=sqrt(RastSQR)
 if (rz>rmax(N1,N2,i1,i2)) then
  rz2=1.0/RastSQR
  rz6=rz2*rz2*rz2
  rz7=rz6/rz
  expp=exp(alf(N1,N2,i1,i2)*(1.0-rz/siga(N1,N2,i1,i2)))
  LJatom=bk1(N1,N2,i1,i2)*expp+bk2(N1,N2,i1,i2)*rz6
  VirLJatom=(dbk1(N1,N2,i1,i2)*rz7+dbk2(N1,N2,i1,i2)*expp)*rz
 else
  LJatom=999999999999.0
  VirLJatom=999999999999.0
 endif
endif
if (ptip==5) then
 if ((RastSQR>p5_l2).and.(RastSQR<p5_dl2)) then
  LJatom=p5_e
  VirLJatom=0.0
  !print *,RastSQR,p5_l2
  !pause
 else
  rz2=siga(N1,N2,i1,i2)*siga(N1,N2,i1,i2)/RastSQR
  rz6=rz2*rz2*rz2
  LJatom=4.0*epsa(N1,N2,i1,i2)*(rz6*rz6-rz6)
  VirLJatom=epsa(N1,N2,i1,i2)*(-rz6*rz6 + 0.5*rz6)*48.0
 endif
endif
if (ptip==7) then
 rz=sqrt(RastSQR)
 rz6=(1.0/rz)**modif
 expp=exp(alf(N1,N2,i1,i2)*(1.0-rz/siga(N1,N2,i1,i2)))
 !ekt=ekt+(kk1(i,j)*expp+kk2(i,j)*expp+kk3(i,j)+kk4(i,j))*rz6
 LJatom=(kkm1(N1,N2,i1,i2)*expp+kkm2(N1,N2,i1,i2))*rz6
 VirLJatom=dkkm1(N1,N2,i1,i2)*(modif/(modif+alf(N1,N2,i1,i2))*expp-1.0)*rz6+&
 &dkkm2(N1,N2,i1,i2)*rz6*expp*rz
 !print *,rz,rz6,LJatom,VirLJatom
!pause
endif
if (ptip==9) then
 rz=sqrt(RastSQR)
 expp=(siga(N1,N2,i1,i2)/rz)**alf(N1,N2,i1,i2)    !не экспонента
 expp6=(siga(N1,N2,i1,i2)/rz)**6
 LJatom=epsa(N1,N2,i1,i2)*dopeps(N1,N2,i1,i2)*(expp-expp6)
 VirLJatom=epsa(N1,N2,i1,i2)*dopeps(N1,N2,i1,i2)*(expp6*6.0 - alf(N1,N2,i1,i2)*expp)

endif
end subroutine
!------------------------------------------------------------------------------------------------
!************************************************************************************************
!************************************************************************************************
!*****************************************************
!     FAST FOURIER ROUTINE - IN ONE DIMENSION        *
!*****************************************************
!
      SUBROUTINE FFS (KAM,A,B,AA,BB,NMAIN,delta)
!  A(K) - original functon
!  B(K) - it's image
!  B(K) = COEFF * DR * SUMMA ( SIN(I*K*PI/N) * A(I) )
!  N = 2**NM
!  DK*DR=PI/N
!  COEFF = 4 * PI    ... DIRECT  TRANS. KAM > 0
!        = 2 / PI**2 ... INVERS. TRANS. KAM < 0
!          !!! A(N) MUST BE ZERO !!!
!
integer NMAIN,nm,n, KAM
real(8) PI,delta
parameter ( PI=3.141592653589793)
integer(4) i,j,L,k,L1,L2
real(8) A(1),B(1),AA(1),BB(1)
real(8) P1,P2,P3,P4,P5,P6,P7,P8,COEF

NM=NMAIN
N=2**NMAIN
AA(1)=0.0
AA(N+1)=0.0
BB(1)=0.0
BB(N+1)=0.0
DO I=2,N
 AA(I)=A(I-1)
 AA(N+I)=0.0
 BB(I)=0.0
 BB(N+I)=0.0
enddo
NM=NM+1
N=N+N
J=1
L=N-1
DO I=1,L
 IF (I>=J) then !GOTO 20
  K=N/2
  else
  P1=AA(J)
  AA(J)=AA(I)
  AA(I)=P1
  K=N/2
 endif
 do while (K<J)
  J=J-K
  K=K/2
 enddo
 J=J+K
enddo
DO L=1,NM
 L1=2**L
 L2=L1/2
 P1=1.0
 P2=0.0
 P3=COS(PI/L2)
 P4=-SIN(PI/L2)
 DO  J=1,L2
  DO  I=J,N,L1
   P7=AA(I+L2)
   P8=BB(I+L2)
   P5=P7*P1-P8*P2
   P6=P7*P2+P8*P1
   AA(I+L2)=AA(I)-P5
   AA(I)=AA(I)+P5
   BB(I+L2)=BB(I)-P6
   BB(I)=BB(I)+P6
  enddo
  P5=P1
  P1=P1*P3-P2*P4
  P2=P5*P4+P2*P3
 enddo
enddo
M=NM-1
N=N/2

IF  (KAM==+1) COEF=-4.0*PI*Delta
IF  (KAM==-1) COEF=-0.5/PI/PI*Delta
DO I=2,N
 B(I-1)=COEF*BB(I)
enddo
B(N)=0.0
end subroutine
subroutine calc_cav()
use dannie
!do i=1,ckol,1 !полная корелляционная функция
! do j=i,ckol,1
!  write (rdffile,'(5a)') 'H-',trim(adjustl(labckol(i))),'-',trim(adjustl(labckol(j))),'.txt'
!  open(23, file=rdffile)
!  do k=1,HistRaz,1
!   total_c_f((i-1)*ckol+j,k)=VirGrA((i-1)*ckol+j,k)-1.0
!   write(23,'(f15.10,a1,f15.10)') Grx(k),';', total_c_f((i-1)*ckol+j,k)
!  enddo
!  close(23)
! enddo
!enddo
deltak=1.0/deltar*ChPI/2**stepen
do i=1,ckol,1 !прямая корелляционная функция
 do j=i,ckol,1
  do k=1,HistRaz,1
   AFFT(k)=total_c_f((i-1)*ckol+j,k)*Grx(k)
  enddo
  !AFFT(HistRaz)=0.0
  call FFS(1,AFFT,BFFT,WorkA,WorkB,stepen,deltar)
  do k=1,HistRaz,1
   AFFT(k)=BFFT(k)/(1.0+RoNV*BFFT(k)/(Grx(k)/deltar*deltak))!k/deltak)
  enddo
  call FFS(-1,AFFT,BFFT,WorkA,WorkB,stepen,deltak)
  do k=1,HistRaz,1
   direct_c_f((i-1)*ckol+j,k)=BFFT(k)/Grx(k)
  enddo
 enddo
enddo
do i=1,ckol,1 !вывод прямая корелляционная функция
 do j=i,ckol,1
  write (rdffile,'(5a)') 'C-',trim(adjustl(labckol(i))),'-',trim(adjustl(labckol(j))),'.txt'
  open(23, file=rdffile)
  do k=1,HistRaz,1
   write(23,'(f15.10,a1,f15.10)') Grx(k),';', direct_c_f((i-1)*ckol+j,k)
  enddo
  close(23)
 enddo
enddo
end subroutine

subroutine calc_bridge !,\бридж функционал и непрямая корреляционная
use dannie
do i=1,ckol,1
 do j=1,ckol,1
  write (rdffile,'(5a)') 'Br-',trim(adjustl(labckol(i))),'-',trim(adjustl(labckol(j))),'.txt'
  open(25,file=rdffile)
  do k=1,HistRaz,1
   bridge_f((i-1)*ckol+j,k)=direct_c_f((i-1)*ckol+j,k)+cavity_out((i-1)*ckol+j,k)-&
   &total_c_f((i-1)*ckol+j,k)
   write(25,'(f15.10,a1,f15.10)') Grx(k),';', bridge_f((i-1)*ckol+j,k)
  enddo
  close(25)
  write (rdffile,'(5a)') 'Gamma-',trim(adjustl(labckol(i))),'-',trim(adjustl(labckol(j))),'.txt'
  open(25,file=rdffile)
  do k=1,HistRaz,1
   indirect_c_f((i-1)*ckol+j,k)=total_c_f((i-1)*ckol+j,k)-direct_c_f((i-1)*ckol+j,k)
   write(25,'(f15.10,a1,f15.10)') Grx(k),';', indirect_c_f((i-1)*ckol+j,k)
  enddo
  close(25)
 enddo
enddo
end subroutine

subroutine Saveme()
use dannie
character(2200) f_input
character(2200) f_xa
character(2200) f_x
character(2200) f_summs
character(2200) f_ntip
namelist /main_save/RoNV,Nv,sproc,storona,Naz,tipResh,Temp,RadCut,TipCut, TipLat,Neqv,&
&Ntot,Nnew,Nchan,drtail,Nstep,drtail,workdir,ena,npr,ptip,calc_rdf,calc_hend

write(f_input,'(a,a,a,a1,i2,a1,f8.4,a1,f8.4,a)') trim(adjustl(workdir)),'/save/',&
&NPR,storona,'-',RoNV,'-',Temp, 'input.txt'
write(f_x,'(a,a,a1,i2,a1,f8.4,a1,f8.4,a)') trim(adjustl(workdir)),'/save/',NPR,storona&
&,'-',RoNV,'-',Temp, 'mol.txt'
write(f_xa,'(a,a,a1,i2,a1,f8.4,a1,f8.4,a)') trim(adjustl(workdir)),'/save/',NPR,storona,&
&'-',RoNV,'-',Temp, 'atom.txt'
write(f_summs,'(a,a1,i2,a1,f8.4,a1,f8.4,a)') trim(adjustl(workdir)),'/save/',NPR,storona&
&,'-',RoNV,'-',Temp, 'sums.txt'
write(f_ntip,'(a,a,a1,i2,a1,f8.4,a1,f8.4,a)') trim(adjustl(workdir)),'/save/',NPR,storona,&
&'-',RoNV,'-',Temp, 'ntip.txt'
open(28,file=f_ntip)
open(29,file=f_xa)
open(30,file=f_x)
do i=1,N,1
 write(28,'(i3)') Ntip(i)
 write(30,'(f30.20,f30.20,f30.20)') x(i),y(i),z(i)
 do j=1,Nm(i),1
  write(29,'(f30.20,f30.20,f30.20)') xa(i,j),ya(i,j),za(i,j)
 enddo
enddo

open (31,file=f_input)
 write(31,main_save)
close(31)
end subroutine

subroutine smooth_rdf()
use dannie
integer(4) k_smooth,s_num
real(8) s_sum
k_smooth=10
do inom=1,ckol,1
 do jnom=inom,ckol,1
  do j=1,k_smooth,1
   do i=1,HistRaz,1
    s_sum=0.0
    s_num=0
    !граничные
    if (i>HistRaz-ceiling(rGrcut-0.5*maxsig)) then
     do s_i=0,ceiling(rGrcut-0.5*maxsig)-1,1
      s_sum=s_sum+total_c_f((inom-1)*ckol+jnom,i-s_i)
      s_num=s_num+1
     enddo
     smooth(i)=s_sum/float(s_num)
    else if (deltar*float(i)>1.0*maxsig) then !неграничные
     do s_i=-(ceiling(deltar*(float(i))/maxsig)-1),ceiling(deltar*(float(i))/maxsig)-1,1
      s_sum=s_sum+total_c_f((inom-1)*ckol+jnom,i+s_i)
     enddo
     smooth(i)=s_sum/float(2*(ceiling(deltar*(float(i))/maxsig)-1)+1)
    else
     smooth(i)=total_c_f((inom-1)*ckol+jnom,i)
    endif
   enddo
   do i=1,HistRaz,1
    total_c_f((inom-1)*ckol+jnom,i)=smooth(i)
   enddo
  enddo
 enddo
enddo

do i=1,ckol,1 !полная корелляционная функция
 do j=i,ckol,1
  write (rdffile,'(5a)') 'HS-',trim(adjustl(labckol(i))),'-',trim(adjustl(labckol(j))),'.txt'
  open(23, file=rdffile)
  do k=1,HistRaz,1
   write(23,'(f15.10,a1,f15.10)') Grx(k),';', total_c_f((i-1)*ckol+j,k)
  enddo
  close(23)
 enddo
enddo

do inom=1,ckol,1
 do jnom=inom,ckol,1
  write (rdffile,'(a4,a,a1,a,a4)') 'RDFS-',trim(adjustl(labckol(inom))),'-',trim(adjustl(labckol(jnom))),'.txt'
  open(23, file=rdffile)
  do i=1,HistRaz,1
   write(23,'(f15.10,a1,f15.10)') Grx(i),';', VirGrA((inom-1)*ckol+jnom,i)
  enddo
  close(23)
 enddo
enddo
end subroutine

subroutine outfile(out_int)
use dannie
integer(4) out_int

open(13, file='OUT.txt')
if (CALC_RDF==1) then
 write (13,'(a)') 'Histogram method (TotalCF, HalfCavityCF, DirectCF)'
 if (CALC_HEND==1) then
  write (13,'(a)') 'Henderson method (FullCavityCF, Chemical Potential)'
 endif
endif
do i=1,Nv,1
 if ((ptip==1).or.(ptip==4).or.(ptip==5).or.(ptip==7)) then
  write (13,'(a,a,a3,i3)') 'centers in ',trim(adjustl(MMName(i))),' ' ,Nm(i)
  write (13,'(a,a,a3,f7.3,a,i4,a)') 'Mol. fraction of ',trim(adjustl(MMName(i))),': ',&
  & Sproc(i)*100.0,' %  ',kolNv(i), ' particles'
  do j=1,Nm(Nv),1
   write (13,'(a,a,f6.3,a,f7.4,a,a,f10.5,a)') trim(adjustl(labelm(i,j))), ' α ', alfa(i,j),&
   &' rm ', sigma(i,j), trim(razm_rast), ' ε/kb ', epsi(i,j), trim(razm_tem)
  enddo
 endif
 if ((ptip==2).or.(ptip==3)) then
  write (13,'(a,a,a3,i3)') 'centers in ',trim(adjustl(MMName(i))),' ' ,Nm(i)
  write (13,'(a,a,a3,f7.3,a,i4,a)') 'Mol. fraction of ',trim(adjustl(MMName(i))),': ',&
  & Sproc(i)*100.0,' %  ',kolNv(i), ' particles'
  do j=1,Nm(Nv),1
   write (13,'(a,a,f7.4,a,a,f10.5,a)') trim(adjustl(labelm(i,j))),&
   &' σ ', sigma(i,j), trim(razm_rast), ' ε/kb ', epsi(i,j), trim(razm_tem)
  enddo
 endif
 if (ptip==9) then
  write (13,'(a,a,a3,i3)') 'centers in ',trim(adjustl(MMName(i))),' ' ,Nm(i)
  write (13,'(a,a,a3,f7.3,a,i4,a)') 'Mol. fraction of ',trim(adjustl(MMName(i))),': ',&
  & Sproc(i)*100.0,' %  ',kolNv(i), ' particles'
  do j=1,Nm(Nv),1
   write (13,'(a,a,f6.3,a,f7.4,a,a,f10.5,a)') trim(adjustl(labelm(i,j))), ' n ', alfa(i,j),&
   &' σ ', sigma(i,j), trim(razm_rast), ' ε/kb ', epsi(i,j), trim(razm_tem)
  enddo
 endif
 !
if (ptip==1) then
 write (13,'(a)') trim(adjustl('Karr-Konowalow'))
endif
if (ptip==2) then
 write (13,'(a)') trim(adjustl('Lennard-Jones'))
endif
if (ptip==3) then
 write (13,'(a)') trim(adjustl('Lennard-Jones+Potential Well'))
 write (13,'(a,f15.10)') trim(adjustl('l*=')),lp
 write (13,'(a,f15.10)') trim(adjustl('De=')),Dep
endif
if (ptip==4) then
 write (6,'(a)') trim(adjustl('Potential:   Buckingham (exp-6)'))
endif
if (ptip==5) then
 write (13,'(a)') 'Potential: Karr-Konowalow + square-well'
 write (13,'(a,f15.7)') 'e*: ', p5_e
 write (13,'(a,f15.7)') 'l:  ', p5_l
 write (13,'(a,f15.7)') 'l+Δl', p5_l+p5_dl
endif
if (ptip==7) then
 write (13,'(a)') 'Modifed Karr-Konowalow'
 write (13, '(a,f6.4)') 'Attractive parameter: ',modif
endif
if (ptip==9) then
 write (13,'(a)') 'n-6 Lennard-Jones'
endif
!
 do j=1,ckol,1
  write (13,'(a,a3,i5)') trim(adjustl(labckol(j))),' - ', Ncrdf(i,j)
 enddo
 !длинны связи
 if (bond_kol(i)>0) then
  write(13,'(a)') 'Start bond length: '
  do j=1, bond_kol(i),1
   x1valent=xm(i,bond_n(i,j,1))-xm(i,bond_n(i,j,2))
   y1valent=ym(i,bond_n(i,j,1))-ym(i,bond_n(i,j,2))
   z1valent=zm(i,bond_n(i,j,1))-zm(i,bond_n(i,j,2))
   write(13,'(a,a,a,a,f10.5,a)') trim(adjustl(labelm(i,bond_n(i,j,1)))), '-',trim(adjustl(labelm&
   &(i,bond_n(i,j,2)))), '  : ', sqrt(x1valent*x1valent+&
   &y1valent*y1valent+z1valent*z1valent), razm_rast
  enddo
 else
  write (13,'(a)') 'Chemical bonds: None '
 endif
 !валентные углы
 if (bind_val(i)>0) then
 write(13,'(a)') 'Start valent angles: '
 do j=1,bind_val(i),1
  x1valent=xm(i,bind_vn(i,j,1))-xm(i,bind_vn(i,j,2))
  x2valent=xm(i,bind_vn(i,j,3))-xm(i,bind_vn(i,j,2))
  y1valent=ym(i,bind_vn(i,j,1))-ym(i,bind_vn(i,j,2))
  y2valent=ym(i,bind_vn(i,j,3))-ym(i,bind_vn(i,j,2))
  z1valent=zm(i,bind_vn(i,j,1))-zm(i,bind_vn(i,j,2))
  z2valent=zm(i,bind_vn(i,j,3))-zm(i,bind_vn(i,j,2))
  valcos=(x1valent*x2valent+y1valent*y2valent+z1valent*z2valent)/&
  &sqrt(x1valent*x1valent+y1valent*y1valent+z1valent*z1valent)/&
  &sqrt(x2valent*x2valent+y2valent*y2valent+z2valent*z2valent)
  !print *, sqrt(x1valent*x1valent+y1valent*y1valent+z1valent*z1valent)
  !print *, sqrt(x2valent*x2valent+y2valent*y2valent+z2valent*z2valent)
  !расчет валентных углов
  write(13,'(6a,f12.7,a)') trim(adjustl(labelm(i,bind_vn(i,j,1)))), '-',&
  trim(adjustl(labelm(i,bind_vn(i,j,2)))), '-',  &
  &trim(adjustl(labelm(i,bind_vn(i,j,3)))), ' :',  acos(valcos)/ChPI*180.0, ' º'
 enddo
 else
  write (13,'(a)') 'Valent angles: None '
 endif
 if (bind_kol(i)>0) then
  write(13, '(a)') 'Valent angle not rigid'
  do j=1,bind_kol(i),1
   write(13,'(6a)') trim(adjustl(labelm(i,bind1(i,j)))), '-',trim(adjustl(labelm(&
   &i,bind_cent(i,j)))),'-',trim(adjustl(labelm(i,bind2(i,j)))), ' :'
   write(13, '(a,f20.10)') '                       k : ', bind_k(i,j)*2.0
   write(13, '(a,f20.10)') '       eqvilibrium angle : ', bind_fi(i,j)/ChPI*180.0
  enddo
  !ЗАряды
 endif
 if (charg_kol(i)>0) then
  write(13,'(a)') 'Charges :'
  if (charg_kol(i)>0) then
   do j=1,charg_kol(i),1
    write(13,'(f15.10)') ch_q(i,j)
   enddo
  endif
 else
  write (13,'(a)') 'Charges: None'
 endif
!
enddo
write (13,'(a,i5)') trim(adjustl('molecules in volume: ')), N
write (13,'(a,f10.4,a)') trim(adjustl('Mol. density: ')),RoNVML, trim(razm_ro)
write (13,'(a,f10.4,a2)') trim(adjustl('Temperature: ')), TempK, trim(razm_tem)
write (13,'(a,i10)') trim(adjustl('equilibrium steps: ')), Neqv
write (13,'(a,i10)') trim(adjustl('steps per sample: ')), Nnew
write (13,'(a,i10)') trim(adjustl('total steps: ')), Ntot
if (TipLat==1) then
 if (tipResh==1) then
 write (13,'(a)') trim(adjustl('Start lattice - body-centered lattice'))
 else if (tipResh==2) then
 write (13,'(a)') trim(adjustl('Start latice - face-centered lattice'))
 end if
end if
write (13,'(a,f15.8,a)') trim(adjustl('Cut radius: ')), RadCut, trim(razm_rast)
!
write (13, '(a)') '------Energy-----'
write (13, '(a20,f20.10,a)') 'Kinetic Energy ', enid, razm_en
write (13, '(a20,f20.10,a)') 'Tail Energy ', etail1, razm_en
write (13, '(a20,f20.10,a)') 'Potential Energy ', outsenout, razm_en
write (13, '(a20,f20.10,a)') 'Triple Energy ', tr_out, razm_en
write (13, '(a20,f20.10,a)') 'Total Energy ', outsenout+enid, razm_en
write (13, '(a20,f20.10,a)') 'Energy Difference ', o2enp, razm_en
write (13, '(a20,f30.15,a)') 'Average Heat Capacity: ', OutCvOut, razm_cv
write (13, '(a20,f30.15,a)') 'Heat Capacity Difference', o2cv, razm_cv
write (13, '(a)') '------Pressure-----'
write (13, '(a20,f20.10,a)') 'Pressure ', outdavlout, razm_davl
write (13, '(a20,f20.10,a)') 'Triple Pressure ', tr_dout, razm_davl
write (13, '(a20,f20.10,a)') 'Tail Pressure', ptail1, razm_davl
write (13, '(a20,f20.10,a)') 'Pressure Difference ', o2davl, razm_davl
if (ena==1) then
 write(13,'(a20,f20.10,a)') 'Zv ', outdavlout/RoNVML/TempK
else
 write (13, '(a20,f20.10,a)') 'Zv ', outdavlout*100.0/RoNVML/8.3145107/TempK
endif

write (13,'(a,f20.10,a,f20.10)') 'isotermal compressibility: ', bette, ' bT-1  ', betminus1
write (13, '(a)') '------Chemical Potential-----'
do i=1,Nv,1
 write (13, '(a,a2,f20.10,a)')trim(adjustl(MMName(i))),'- ', mu2_out(i), razm_mu
enddo
do i=1,ckol,1
  do j=i,ckol,1
   write(13,'(a,a,a,a,a,f20.10)') 'Coordination number ',trim(adjustl(labckol(i))),'-',trim(adjustl&
   &(labckol(j))),' :  ', kc(i,j)
  enddo
 enddo
write (13,'(a,i10,a,i10)') 'accepted ', assept, ' not accepted ', notassept
!write (13, '(a,f15.7,a)') 'computation time ', finish-start ,'sec'
!if (out_int==2) then
! write (13, '(a,f15.7,a)') 'computation time ', finish-start ,'sec'
! hour_t=int(((finish-start)-mod(finish-start,3600.0))/3200.0)
! min_t=int(((finish-start)-Hour_t*3200.0-mod(finish-start,60.0))/60.0)
! sec_t=mod(finish-start,60.0)
! write (13, '(a,i3,a,i2,a,f15.7,a)') 'computation time ',&
! & hour_t ,' hour ', min_t, ' minute ', sec_t ,'sec'
! write (13, '(a)') ' calculation is finished '
! !write (13,'(5f30.15)') KonSt, finish-start, sen/Meqvil,sen/Meqvil+enid, sdavl/Meqvil
!else
 write (13, '(a,f15.7,a)') 'computation time ', time1-start ,'sec'
 sec_t=mod(time1-start,60.0)
 min_t=int(mod(((time1-start)-sec_t),3600.0)/60)
 hour_t=int(((time1-start)-min_t-sec_t)/3600)
 !hour_t=int(((time1-start)-mod(time1-start,3600.0))/3200.0)
 !min_t=int(((time1-start)-Hour_t*3200.0-mod(time1-start,60.0))/60.0)

 write (13, '(a,i3,a,i2,a,f15.7,a)') 'computation time ',&
 & hour_t ,' hour ', min_t, ' minute ', sec_t ,'sec'
! write (13, '(a,f15.7,a)') 'computation time for sample ', time1-time2 ,'sec'
!endif
close(13)
!
open(13, file='OUTLATA.xyz') !вывод относительного положения атомов
write(13,'(i3)') totalNa
write(13,'(a52)') 'XYZ file generated by MCP : coordinates in Angstrom'
do i=1,N,1
 do j=1,Na(i),1
  write (13,'(a2,3f20.10)') labela(i,j),(x(i)+xa(i,j))*bs,(y(i)+ya(i,j))*bs,(z(i)+za(i,j))*bs
 enddo
enddo
close(13)
!
end subroutine

subroutine bind_ch(Nbind)
use dannie
integer(4) Nbind,ib,jto
real(8) bind_rnd
real(8) sin1,sin2,sin3,sin4,cos1,cos2,cos3,cos4
do ib=1,bind_kol(Ntip(Nbind)),1
  !ставим координаты молекулы в 000
 xzero=xa(Nbind,bind_cent(Ntip(Nbind),ib))  !определение
 yzero=ya(Nbind,bind_cent(Ntip(Nbind),ib))
 zzero=za(Nbind,bind_cent(Ntip(Nbind),ib))
 do i=1,Na(Nbind),1             !перенос молекулы в центр
  xa(Nbind,i)=xa(Nbind,i)-xzero
  ya(Nbind,i)=ya(Nbind,i)-yzero
  za(Nbind,i)=za(Nbind,i)-zzero
 enddo
 ! прямая пеерестановка
 !
!   xv1=xa(Nbind,1)-xa(Nbind,2)
!  yv1=ya(Nbind,1)-ya(Nbind,2)
!  zv1=za(Nbind,1)-za(Nbind,2)
!  xv2=xa(Nbind,3)-xa(Nbind,2)
!  yv2=ya(Nbind,3)-ya(Nbind,2)
!  zv2=za(Nbind,3)-za(Nbind,2)
!  !print *, xv1,xv2,yv1,yv2,zv1,zv2
!  r1_bind=sqrt(xv1*xv1+yv1*yv1+zv1*zv1)
!  r2_bind=sqrt(xv2*xv2+yv2*yv2+zv2*zv2)
!  print *, (xv1*xv2+yv1*yv2+zv1*zv2)/(r1_bind*r2_bind)
 !
 !после этого центр ориентирован на оси y (z=0)
 if (ya(Nbind,bind1k(Ntip(Nbind),ib,1))==0.0) then  !на случай нулевых значений
  sin1=0.0
  cos1=1.0
 else
  sin1=ya(Nbind,bind1k(Ntip(Nbind),ib,1))/sqrt(ya(Nbind,bind1k(Ntip(Nbind),ib,1))&
  &*ya(Nbind,bind1k(Ntip(Nbind),ib,1))+za(Nbind,bind1k(Ntip(Nbind),ib,1))*za(Nbind&
  &,bind1k(Ntip(Nbind),ib,1)))
  cos1=za(Nbind,bind1k(Ntip(Nbind),ib,1))/sqrt(ya(Nbind,bind1k(Ntip(Nbind),ib,1))&
  &*ya(Nbind,bind1k(Ntip(Nbind),ib,1))+za(Nbind,bind1k(Ntip(Nbind),ib,1))*za(Nbind&
  &,bind1k(Ntip(Nbind),ib,1)))
 endif
 do jto=1,Na(Nbind),1
  yv=-za(Nbind,jto)*sin1+ya(Nbind,jto)*cos1
  zv=za(Nbind,jto)*cos1+ya(Nbind,jto)*sin1
  ya(Nbind,jto)=yv
  za(Nbind,jto)=zv
 enddo
 do jto=1,charg_kol(Ntip(Nbind)),1              !перенос зарядов
  yv=-ch_az(Nbind,jto)*sin1+ch_ay(Nbind,jto)*cos1
  zv=ch_az(Nbind,jto)*cos1+ch_ay(Nbind,jto)*sin1
  ch_ay(Nbind,jto)=yv
  ch_az(Nbind,jto)=zv
 enddo
 !после этого центр ориентирован на оси x (y=0,z=0)
 if (za(Nbind,bind1k(Ntip(Nbind),ib,1))==0.0) then
  sin2=0.0
  cos2=1.0
 else
  sin2=za(Nbind,bind1k(Ntip(Nbind),ib,1))/sqrt(xa(Nbind,bind1k(Ntip(Nbind),ib,1))&
  &*xa(Nbind,bind1k(Ntip(Nbind),ib,1))&
  &+za(Nbind,bind1k(Ntip(Nbind),ib,1))*za(Nbind,bind1k(Ntip(Nbind),ib,1)))
  cos2=xa(Nbind,bind1k(Ntip(Nbind),ib,1))/sqrt(xa(Nbind,bind1k(Ntip(Nbind),ib,1))*&
  &xa(Nbind,bind1k(Ntip(Nbind),ib,1))&
  &+za(Nbind,bind1k(Ntip(Nbind),ib,1))*za(Nbind,bind1k(Ntip(Nbind),ib,1)))
 endif
 do jto=1,Na(Nbind),1
  xv=xa(Nbind,jto)*cos2+za(Nbind,jto)*sin2
  zv=-xa(Nbind,jto)*sin2+za(Nbind,jto)*cos2
  xa(Nbind,jto)=xv
  za(Nbind,jto)=zv
 enddo
 do jto=1,charg_kol(Ntip(Nbind)),1              !перенос зарядов
  xv=ch_ax(Nbind,jto)*cos2+ch_az(Nbind,jto)*sin2
  zv=-ch_ax(Nbind,jto)*sin2+ch_az(Nbind,jto)*cos2
  ch_ax(Nbind,jto)=xv
  ch_az(Nbind,jto)=zv
 enddo
!передвигаем другую сторону
! после этого (z=0)
 if (za(Nbind,bind2k(Ntip(Nbind),ib,1))==0.0) then
  sin3=0.0
  cos3=1.0
 else
  sin3=za(Nbind,bind2k(Ntip(Nbind),ib,1))/sqrt(ya(Nbind,bind2k(Ntip(Nbind),ib,1))&
  &*ya(Nbind,bind2k(Ntip(Nbind),ib,1))&
  &+za(Nbind,bind2k(Ntip(Nbind),ib,1))*za(Nbind,bind2k(Ntip(Nbind),ib,1)))
  cos3=ya(Nbind,bind2k(Ntip(Nbind),ib,1))/sqrt(ya(Nbind,bind2k(Ntip(Nbind),ib,1))&
  &*ya(Nbind,bind2k(Ntip(Nbind),ib,1))&
  &+za(Nbind,bind2k(Ntip(Nbind),ib,1))*za(Nbind,bind2k(Ntip(Nbind),ib,1)))
 endif
 do jto=1,Na(Nbind),1
  yv=ya(Nbind,jto)*cos3+za(Nbind,jto)*sin3
  zv=-ya(Nbind,jto)*sin3+za(Nbind,jto)*cos3
  ya(Nbind,jto)=yv
  za(Nbind,jto)=zv
 enddo
 do jto=1,charg_kol(Ntip(Nbind)),1
  yv=ch_ay(Nbind,jto)*cos3+ch_az(Nbind,jto)*sin3
  zv=-ch_ay(Nbind,jto)*sin3+ch_az(Nbind,jto)*cos3
  ch_ay(Nbind,jto)=yv
  ch_az(Nbind,jto)=zv
 enddo
! крутим часть молекулы относительно оси z
 !bind_rnd=
 !bind_u(Nbind,ib)=bind_u(Nbind,ib)+bind_rnd !измеряем потом
 !print *,bind_rnd
 sin4=BindDeg*(getrand()-0.5)*2.0
 cos4=sqrt(1.0-sin4*sin4)
 do jto=1,bind1(Ntip(Nbind),ib),1
  yv=ya(Nbind,bind1k(Ntip(Nbind),ib,jto))*cos4+xa(Nbind,bind1k(Ntip(Nbind),ib,jto))*sin4
  xv=-ya(Nbind,bind1k(Ntip(Nbind),ib,jto))*sin4+xa(Nbind,bind1k(Ntip(Nbind),ib,jto))*cos4
  xa(Nbind,bind1k(Ntip(Nbind),ib,jto))=xv
  ya(Nbind,bind1k(Ntip(Nbind),ib,jto))=yv
 enddo
 !поворачиваем часть молекулы

 !крутим обратно
 do jto=1,Na(Nbind),1
  yv=ya(Nbind,jto)*cos3-za(Nbind,jto)*sin3
  zv=ya(Nbind,jto)*sin3+za(Nbind,jto)*cos3
  za(Nbind,jto)=zv
  ya(Nbind,jto)=yv
 enddo
 do jto=1,charg_kol(Ntip(Nbind)),1          !поворачиваем заряд
  yv=ch_ay(Nbind,jto)*cos3-ch_az(Nbind,jto)*cos3
  zv=ch_ay(Nbind,jto)*sin3+ch_az(Nbind,jto)*cos3
  ch_az(Nbind,jto)=zv
  ch_ay(Nbind,jto)=yv
 enddo
 !крутим еще обратно
 do jto=1,Na(Nbind),1
  xv=xa(Nbind,jto)*cos2-za(Nbind,jto)*sin2
  zv=xa(Nbind,jto)*sin2+za(Nbind,jto)*cos2
  za(Nbind,jto)=zv
  xa(Nbind,jto)=xv
 enddo
 do jto=1,charg_kol(Ntip(Nbind)),1          !поворачиваем заряд
  xv=ch_ax(Nbind,jto)*cos2-ch_az(Nbind,jto)*sin2
  zv=ch_ax(Nbind,jto)*sin2+ch_az(Nbind,jto)*cos2
  ch_az(Nbind,jto)=zv
  ch_ax(Nbind,jto)=xv
 enddo
 do jto=1,Na(Nbind),1
  yv=za(Nbind,jto)*sin1+ya(Nbind,jto)*cos1
  zv=za(Nbind,jto)*cos1-ya(Nbind,jto)*sin1
  ya(Nbind,jto)=yv
  za(Nbind,jto)=zv
 enddo
 do jto=1,charg_kol(Ntip(Nbind)),1
  yv=ch_az(Nbind,jto)*sin1+ch_ay(Nbind,jto)*cos1
  zv=ch_az(Nbind,jto)*cos1-ch_ay(Nbind,jto)*sin1
  ch_ay(Nbind,jto)=yv
  ch_az(Nbind,jto)=zv
 enddo
 do i=1,Na(Nbind),1
  xa(Nbind,i)=xa(Nbind,i)+xzero
  ya(Nbind,i)=ya(Nbind,i)+yzero
  za(Nbind,i)=za(Nbind,i)+zzero
 enddo
 do i=1,charg_kol(Ntip(Nbind)),1
  ch_ax(Nbind,i)=ch_ax(Nbind,i)+xzero
  ch_ay(Nbind,i)=ch_ay(Nbind,i)+yzero
  ch_az(Nbind,i)=ch_az(Nbind,i)+zzero
 enddo
  !
!  xv1=xa(Nbind,1)-xa(Nbind,2)
!  yv1=ya(Nbind,1)-ya(Nbind,2)
!  zv1=za(Nbind,1)-za(Nbind,2)
!  xv2=xa(Nbind,3)-xa(Nbind,2)
!  yv2=ya(Nbind,3)-ya(Nbind,2)
!  zv2=za(Nbind,3)-za(Nbind,2)
!  !print *, xv1,xv2,yv1,yv2,zv1,zv2
!  r1_bind=sqrt(xv1*xv1+yv1*yv1+zv1*zv1)
!  r2_bind=sqrt(xv2*xv2+yv2*yv2+zv2*zv2)
!  print *, (xv1*xv2+yv1*yv2+zv1*zv2)/(r1_bind*r2_bind)
!  pause
 !
enddo

end subroutine
!подсчет энергии валентного взаимодействия
subroutine bind_energy(Nbind)
use dannie
integer(4) ib
bind_en=0.0
do ib=1,bind_kol(Ntip(Nbind)),1
!
  xv1=xa(Nbind,bind1k(Ntip(Nbind),ib,1))-xa(Nbind,bind_cent(Ntip(Nbind),ib))
  yv1=ya(Nbind,bind1k(Ntip(Nbind),ib,1))-ya(Nbind,bind_cent(Ntip(Nbind),ib))
  zv1=za(Nbind,bind1k(Ntip(Nbind),ib,1))-za(Nbind,bind_cent(Ntip(Nbind),ib))
  xv2=xa(Nbind,bind2k(Ntip(Nbind),ib,1))-xa(Nbind,bind_cent(Ntip(Nbind),ib))
  yv2=ya(Nbind,bind2k(Ntip(Nbind),ib,1))-ya(Nbind,bind_cent(Ntip(Nbind),ib))
  zv2=za(Nbind,bind2k(Ntip(Nbind),ib,1))-za(Nbind,bind_cent(Ntip(Nbind),ib))
  !print *, xv1,xv2,yv1,yv2,zv1,zv2
  r1_bind=(xv1*xv1+yv1*yv1+zv1*zv1)
  r2_bind=(xv2*xv2+yv2*yv2+zv2*zv2)
  cos_bind=(xv1*xv2+yv1*yv2+zv1*zv2)/sqrt(r1_bind*r2_bind)
  bind_u(Nbind,ib)=acos(cos_bind)
!
 bind_en=bind_en+bind_k(Ntip(Nbind),ib)*(bind_u(Nbind,ib)-bind_fi(Ntip(Nbind),ib))&
 &*(bind_u(Nbind,ib)-bind_fi(Ntip(Nbind),ib))
 !print *, bind_u(Nbind,ib)/ChPI*180.0,bind_fi(Ntip(Nbind),ib)/ChPI*180.0
 !print *,bind_en,bind_u(Nbind,ib)-bind_fi(Ntip(Nbind),ib)
 !pause
enddo
end subroutine

subroutine bind_ditr() !распределение по углам
use dannie
integer(4) h_bin
!print *,'ok'
do i=1,N,1 !рассматриваем все молекулы
 do j=1,bind_kol(Ntip(i)),1 !рассматриваем все углы поворота
  h_bin=ceiling(bind_u(i,j)/ChPI*180.0/deg_sh)
  !print *, h_bin
  hist_bind(ntip(i),j,h_bin)=hist_bind(ntip(i),j,h_bin)+1
  hist_norm(ntip(i),j)=hist_norm(ntip(i),j)+1
 enddo
enddo
!print *,'ok'
!вывод
sum_hist=0
do i=1,Nv,1
 do j=1,bind_kol(i),1
  write(bindfile,'(a,i2)') trim(adjustl(MMname(i))),j
  open(45,file=bindfile)
  do k=1,kol_deg-1,1
   write(45,'(f15.10,a1,f15.10)') hist_u(i,j,k),';',float(hist_bind(i,j,k))/float(hist_norm(i,j))
   !проверка
   sum_hist=sum_hist+hist_bind(i,j,k)
  enddo
  !if (sum_hist==hist_norm(i,j)) then !проверка!
  ! !print *,'ok'
  !endif
  close(45)
 enddo
enddo
end subroutine

subroutine ewald_setup(kappa)
use dannie
real(8) b,rxk,ryk,rzk
integer(4) kx,ky,kz,totk
!
b=1.0/4.0/kappa/kappa
totk=0
do kx=1,kmax,1
 rkx=ChPI2*real(kx)
 do ky=-kmax,kmax,1
  rky=ChPI2*real(ky)
  do kz=-kmax,kmax,1
   rkz=ChPI2*real(kz)
   ksq=kx*kx+ky*ky+kz*kz
   if ((ksq< ksqmax).and.(ksq/=0)) then
   totk=totk+1
   !if (totk>maxk) then stop 'kvec is too small'
   rksq=rkx*rkx+rky*rky+rkz*rkz
   kvec(totk)=ChPI2*exp(-b*rksq)/rksq
   endif
  enddo
 enddo
enddo
write( 6, '(a,i5)' ) 'Number of vector: ',totk
end subroutine

subroutine ewald_r(kappa,Nmol)
use dannie
!Рассчитывается для каждой молекулы
do i1=1,charg_kol(Ntip(Nmol)),1   !цикл по всем зарядам молекулы
 !определяем растояние
 do i2=1,N,1   !цикл по всем молекулам
  if (Nmol/=i2) then   !если молекула другая

  endif
 enddo
enddo
end subroutine

subroutine calc_tors(Nmol,Ntrs)
use dannie
integer(4) Nmol,Ntrs
do k=1,Ntors(Ntip(Nmol))
 !определение торсионого угла
 d0x=xa(Nmol,rota(Ntip(Nmol),Ntrs,1))-xa(Nmol,fsta(Ntip(Nmol),Ntrs))
 d0y=ya(Nmol,rota(Ntip(Nmol),Ntrs,1))-ya(Nmol,fsta(Ntip(Nmol),Ntrs))
 d0z=za(Nmol,rota(Ntip(Nmol),Ntrs,1))-za(Nmol,fsta(Ntip(Nmol),Ntrs))
 !
 d1x=xa(Nmol,seca(Ntip(Nmol),Ntrs))-xa(Nmol,fsta(Ntip(Nmol),Ntrs))
 d1y=ya(Nmol,seca(Ntip(Nmol),Ntrs))-ya(Nmol,fsta(Ntip(Nmol),Ntrs))
 d1z=za(Nmol,seca(Ntip(Nmol),Ntrs))-za(Nmol,fsta(Ntip(Nmol),Ntrs))
 !
 d2x=xa(Nmol,srota(Ntip(Nmol),Ntrs,1))-xa(Nmol,seca(Ntip(Nmol),Ntrs))
 d2y=ya(Nmol,srota(Ntip(Nmol),Ntrs,1))-ya(Nmol,seca(Ntip(Nmol),Ntrs))
 d2z=za(Nmol,srota(Ntip(Nmol),Ntrs,1))-za(Nmol,seca(Ntip(Nmol),Ntrs))
 !
 d0d1x=d0y*d1z-d0z*d1y
 d0d1y=d0z*d1x-d0x*d1z
 d0d1z=d0x*d1y-d0y*d1x
 !
 d1d2x=d1y*d2z-d1z*d2y
 d1d2y=d1z*d2x-d1x*d2z
 d1d2z=d1x*d2y-d1y*d2x
 !
 TorsCosFi=(d0d1x*d1d2x+d0d1y*d1d2y+d0d1z*d1d2z)/sqrt(d0d1x*d0d1x+d0d1y*d0d1y+&
 &d0d1z*d0d1z)/sqrt(d1d2x*d1d2x+d1d2y*d1d2y+d1d2z*d1d2z)
 TorsFi=acos(TorsCosFi)/ChPI*180.0
 !проверка на 0/180
 if (TorsCosFi==1) then
  if (d0d1x*d1d2x+d0d1y*d1d2y+d0d1z*d1d2z>0.0000001) then
   TorsFi=180.0
  endif
 endif
enddo
end subroutine

subroutine tors_dist()
use dannie
integer(4) tHist
!
open(58,file='starttors.txt')
do i=1,N,1
 do j=1,Ntors(Ntip(i))
  call tors_energy(i)
  call calc_tors(i,j)
  write (58,'(i5,i5,f20.10,f20.10,f20.10)') i,j,TorsCosFi,TorsFi,TorsEnergy
 enddo
enddo
close(58)
!
do i=1,N,1
 do j=1,Ntors(Ntip(i)),1
  tHist=ceiling(TorsU(i,j)/0.25)
  TorsDist(Ntip(i),j,tHist)=TorsDist(Ntip(i),j,tHist)+1.0
 enddo
enddo
end subroutine

!subroutine triple_pm(Nmol)
!use dannie


!end subroutine

subroutine triple_pot(TNmol)
use dannie
integer(4) tr_i,tr_j,tr_k
integer(4) tr_a,tr_b,tr_c
integer(4) TNmol
real(8) tri_px1,tri_py1,tri_pz1
real(8) tri_px2,tri_py2,tri_pz2
real(8) tri_px3,tri_py3,tri_pz3
real(8) tri_r1,tri_r2,tri_r3

real(8) tr3_r1,tr3_r2,tr3_r3
real(8) tr5_r1,tr5_r2,tr5_r3
real(8) tr2_r1,tr2_r2,tr2_r3


Uijk=0.0                    !обнуляем потенциал

do tr_a=1,Na(TNmol),1       !проходимся по центрам молекулы
 do tr_i=1,N               !проходимся по другим молекулам
  if (tr_i/=TNmol) then     !кроме молекулы с которой проверяем
   do tr_b=1,Na(tr_i),1    !по всем центрай другой молекулы
    !
    do tr_j=tr_i+1,N       !ну и по третей
    !
    if (tr_j/=TNmol) then
     do tr_c=1,Na(tr_j)
     !вобщем то само определение потенциала
     tri_px1=(x(TNmol)+xa(TNmol,tr_a))-(x(tr_i)+xa(tr_i,tr_b))+trx_l(tr_i) !!tri_x(TNmol,tr_i*tr_b)
     tri_px2=(x(TNmol)+xa(TNmol,tr_a))-(x(tr_j)+xa(tr_j,tr_c))+trx_l(tr_j) !!tri_x(TNmol,tr_j*tr_c)
     tri_px3=(x(tr_i)+xa(tr_i,tr_b))-(x(tr_j)+xa(tr_j,tr_c))+trx_l(tr_i)+trx_l(tr_j)
     !
     tri_py1=(y(TNmol)+ya(TNmol,tr_a))-(y(tr_i)+ya(tr_i,tr_b))+try_l(tr_i) !!tri_y(TNmol,tr_i*tr_b)
     tri_py2=(y(TNmol)+ya(TNmol,tr_a))-(y(tr_j)+ya(tr_j,tr_c))+trx_l(tr_j) !!tri_y(TNmol,tr_j*tr_c)
     tri_py3=(y(tr_i)+ya(tr_i,tr_b))-(y(tr_j)+ya(tr_j,tr_c))+try_l(tr_i)+try_l(tr_j)
     !
     tri_pz1=(z(TNmol)+ya(TNmol,tr_a))-(z(tr_i)+za(tr_i,tr_b))+trz_l(tr_i)  !!tri_z(TNmol,tr_i*tr_b)
     tri_pz2=(z(TNmol)+za(TNmol,tr_a))-(z(tr_j)+za(tr_j,tr_c))+trz_l(tr_j)  !!tri_z(TNmol,tr_j*tr_c)
     tri_pz3=(z(tr_i)+za(tr_i,tr_b))-(z(tr_j)+za(tr_j,tr_c))+trz_l(tr_i)+trz_l(tr_j)
     !
     tri_r1=sqrt(tri_px1*tri_px1+tri_py1*tri_py1+tri_pz1*tri_pz1)
     tri_r2=sqrt(tri_px2*tri_px2+tri_py2*tri_py2+tri_pz2*tri_pz2)
     tri_r3=sqrt(tri_px3*tri_px3+tri_py3*tri_py3+tri_pz3*tri_pz3)
     !
     tr2_r1=tri_r1*tri_r1
     tr2_r2=tri_r2*tri_r2
     tr2_r3=tri_r3*tri_r3

     tr3_r1=tr2_r1*tri_r1
     tr3_r2=tr2_r2*tri_r2
     tr3_r3=tr2_r3*tri_r3

     tr5_r1=tr3_r1*tr2_r1
     tr5_r2=tr3_r2*tr2_r2
     tr5_r3=tr3_r3*tr2_r3
     !print *, KonSt, tr_i,tr_j
     !print *, trx_l(tr_i),trx_l(tr_j)
     !print *,try_l(tr_i),try_l(tr_j)
     !print *,trz_l(tr_i),trz_l(tr_j)
     !print *,tri_r1,tri_r2,tri_r3
     !

     Uijk=Uijk+tr_v(Ntip(TNmol),tr_a,Ntip(tr_i),tr_b,Ntip(tr_j),tr_c)*(1.0/(tr3_r1*tr3_r2*tr3_r3)) !+&
    ! &3.0*(-tr2_r1+tr2_r2+tr2_r3)*(tr2_r1-tr2_r2+tr2_r3)*(tr2_r1+tr2_r2-tr2_r3)/&
    ! &8.0/(tr5_r1*tr5_r2*tr5_r3))
     !if (Uijk<0) then
     ! print *, Uijk, tri_r1, tri_r2,tri_r3
     ! pause
     !endif
     !print *,Uijk,tr_v(Ntip(TNmol),tr_a,Ntip(tr_i),tr_b,Ntip(tr_j),tr_c),tri_r1,tri_r2,tri_r3

      !  if (TNmol==2) then
       ! pause
        !endif
     enddo
    endif
    enddo
    !
   enddo
  endif
 enddo
enddo
     !print *,Uijk,tr_v(Ntip(TNmol),tr_a,Ntip(tr_i),tr_b,Ntip(tr_j),tr_c), TNmol
      !  pause

end subroutine

!subroutine ch_l(N)
!use dannie
!real(8) l_izm
!do i=1,l_kol(Ntip(N))
! do j=1,l_n1(Ntip(N))
!   !определяем случайное изменение длины
!   l_izm=getrand()*0.05*DlShag
!
! enddo
!enddo
!end subroutine

subroutine modif_move(Nmol)
use dannie
real(8) randalfa !случайный поворот относительно оси
real(8) randbeta !случайный поворот относительно оси
real(8) randgamma !случайный поворот относительно оси
real(8) rast_vr
real(8) sinOx,sinOy,sinOz
real(8) cosOx,cosOy,cosOz

randalfa=getrand()
randbeta=getrand()
randgamma=getrand()
!перемещаем молекулы так чтобы вокруг которой была в центре
x(Nmol)=x(Nmol)-x(near_mol(Nmol))
y(Nmol)=y(Nmol)-y(near_mol(Nmol))
z(Nmol)=z(Nmol)-z(near_mol(Nmol))
!находим растояние вокруг которого надо вращать
rast_vr=1.0+(getrand()-0.5)*0.005
!добавляем случайное изменение
x(Nmol)=x(Nmol)*rast_vr
y(Nmol)=y(Nmol)*rast_vr
z(Nmol)=z(Nmol)*rast_vr
!задаем углы поворотов
sinOx=0.4*(randalfa-0.5)*2.0
cosOx=sqrt(1.0-sinOx*sinOx)
sinOy=0.4*(randbeta-0.5)*2.0
cosOy=sqrt(1.0-sinOy*sinOy)
sinOz=0.4*(randgamma-0.5)*2.0
cosOz=sqrt(1.0-sinOz*sinOz)
xv=x(Nmol)
yv=y(Nmol)
zv=z(Nmol)
!поворачиваем
x(Nmol)=cosOx*cosOy*xv+(cosOx*sinOy*sinOz-sinOx*cosOz)*yv+(cosOx*sinOy*cosOz+sinOx*sinOz)*zv
y(Nmol)=sinOx*cosOy*xv+(sinOx*sinOy*sinOz+cosOx*cosOz)*yv+(sinOx*sinOy*cosOz-cosOx*sinOz)*zv
z(Nmol)=-sinOy*xv+cosOy*sinOz*yv+cosOy*cosOz*zv
!вставляем моекулу назад
x(Nmol)=x(Nmol)+x(near_mol(Nmol))
y(Nmol)=y(Nmol)+y(near_mol(Nmol))
z(Nmol)=z(Nmol)+z(near_mol(Nmol))
call proV(Nmol)


end subroutine


subroutine checknear(nmove)
use dannie
integer(4) cnnv,cnnv2,cnnm,cnnm2
real(8) cnrmin
integer(4) ndis,nmove

cnkol=cnkol+1
ndis=0
do cnnm=1,N
    cnrmin=100000000.0
    do cnnm2=1,N
    if (cnnm/=cnnm2) then
        dmx=(-x(cnnm)+x(cnnm2))
        dmx=abs(dmx)
        if (dmx>KonSt/2.0) then
            dmx=dmx-KonSt/2.0
        endif
        dmy=(-y(cnnm)+y(cnnm2))
        dmy=abs(dmy)
        if (dmy>KonSt/2.0) then
            dmy=dmy-KonSt/2.0
        endif
        dmz=(-z(cnnm)+z(cnnm2))
        dmz=abs(dmz)
        if (dmz>KonSt/2.0) then
            dmz=dmz-KonSt/2.0
        endif
        rastmm=sqrt(dmx*dmx+dmy*dmy+dmz*dmz)
        if (rastmm<cnrmin) then
            cnrmin=rastmm
            near_mol(cnnm)=cnnm2
            near_rast(cnnm)=rastmm
        endif
    endif
    enddo
        if (near_rast(cnnm)<glomin) then
            ndismol=ndismol+1
            ndis=ndis+1
        endif
enddo
write(56,'(i10,a,i6,a,i6)') nmove, ';', ndis, ';', N
!    print *,near_rast
!    print *,near_mol
!    pause
!print *,cnkol, ndismol, ndismol/cnkol/N

end subroutine

subroutine xyzanim()
use dannie
integer(4) xyzi,xyzj,xyzk
    write(59,'(i5)') N
    write(59,'(a)') 'comentline'
do xyzi=1,N
    write(59, '(a,f10.6,f10.6,f10.6)') labela(1,1),x(xyzi),y(xyzi),z(xyzi)
enddo

end subroutine

subroutine densdistr()
use dannie
integer(4) densi
integer(4) denshist
integer(4) denshist1,denshist2
!print *, 'begin'
densNo=densNo+1
!print *,KonSt, zdel
do densi=1,N
    denshist=ceiling((z(densi)+KonSt)/zdel*float(densl))
    denshist1=ceiling((z(densi)+KonSt)/zdel*float(densl1))
    denshist2=ceiling((z(densi)+KonSt)/zdel*float(densl2))
    if (denshist>densl) then
        print *, 'memeror', denshist, z(densi), zdel
        !pause
    endif
    Denssum(denshist)=denssum(denshist)+1.0
    Denssum1(denshist1)=denssum1(denshist1)+1.0
    Denssum2(denshist2)=denssum2(denshist2)+1.0
enddo
open (66,file='densdistr.txt')
open (71,file='densdistr1.txt')
open (73,file='densdistr2.txt')
do densi=1,densl
    write(66,'(2f30.15)') (float(densi)-0.5)*zdel/float(densl),&
    & denssum(densi)/float(densNo)/(zdel/float(densl)*KonSt*KonSt)
enddo
do densi=1,densl1
    write(71,'(2f30.15)') (float(densi)-0.5)*zdel/float(densl1),&
    & denssum1(densi)/float(densNo)/(zdel/float(densl1)*KonSt*KonSt)
enddo
do densi=1,densl2
    write(73,'(2f30.15)') (float(densi)-0.5)*zdel/float(densl2),&
    & denssum2(densi)/float(densNo)/(zdel/float(densl2)*KonSt*KonSt)
enddo

close(66)
close(71)
close(73)

!print *, 'end'

end subroutine
