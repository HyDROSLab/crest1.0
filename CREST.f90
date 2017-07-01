Program CREST
	Implicit None
	Character(80):: fileName, msg, msgD
	Character(80), Allocatable::expMark(:)
	Integer(2):: I, Dcur(1:4), J
	Logical:: fExist

	Character(80):: FDbasic, fDrain, fDpet, fDexport
	Integer(1):: RainStyle, PETstyle, EchoLevel, conExport(0:15), ExpFormat
	Integer(2):: DB(1:4)
	Integer(4):: Nmax, Nwarm, outN
	Integer(4), Allocatable:: expR(:), expC(:)
	Real(4):: dT
	Real(4), Allocatable:: expMean(:,:)
	
	filename="Control.txt"
	Inquire(file=Trim(fileName), exist=fExist)
	If (fExist) Then
		Open(102,file=Trim(fileName), form='formatted')
			Read(102,*) fDbasic, msg
			Read(102,*) fDexport, ExpFormat, msg
			Read(102,*) fDrain, RainStyle, msg
			Read(102,*) fDpet, PETstyle, msg
			Read(102,*) msgD, dT, Nmax, Nwarm, EchoLevel, msg
			Call D_yyyymmddhh(DB, msgD)
			Do I = 0, 15
				Read(102,*) conExport(I), msg
			End Do
			Read(102,*) outN, msg
			Allocate(expR(0:outN), expC(0:outN), expMark(0:outN))
			Do I = 1, outN
				Read(102,*) expR(I), expC(I), expMark(I)
			End Do
			Allocate(expMean(1:(14+outN), 1:Nmax))
		Close(102)

		Call CREST_Excute(FDbasic, fDrain, fDpet, fDexport, RainStyle, PETstyle, &
			&DB, dT, Nmax, Nwarm, EchoLevel, conExport, ExpFormat, outN, expMean, expR, expC)

		!Output Results
		Call Nowmsg(msg)
		fileName=Trim(fDexport) // "CREST." // Trim(msg) //".Mean.txt"
		Open(102, file=fileName, form='formatted')
			Write(102, "(A10, $)") "Date" , "Rain", "PET", "actE", "ExcSur", "ExcBas", "Excess", "StoSur", "StoBas", "Storage", "WU", "WL", "WD", "Runoff", "PmaxC"
			Do I = 1, outN
				 Write(102, "(A10, $)") Trim(expMark(I))
			End Do
			Write(102,"(A5)") "EOL"
			
			Write(102,"(A10, $)") "", "mm/h", "mm/h", "mm/h", "mm/h", "mm/h", "mm/h", "mm/h", "mm/h", "mm/h", "%", "%", "%", "mm/h", "m^3/s"
			Do I = 1, outN
				 Write(102,"(A10, $)") "m^3/s"
			End Do
			Write(102,"(A5)") "EOL"
				
			Dcur=DB
			Do I = Nwarm + 1, Nmax
				Call D2yyyymmddhh(Dcur,msg)
				 Write(102,"(A10, $)") Trim(msg)
				 Do J = 1, 14+outN
					  Write(102,"(F10.3, $)") expMean(J, I)
				 End Do
				 Write(102,"(A5)") "EOL"
				 Call Date_Change(Dcur, 4, dT)
			End Do
		Close(102)
	Else
		Write(*,*) "Main File Missing " // fileName
	End If
End Program


! DEM is Mask File
! Unit of DT is Hours£¬Unit of Precipitation and PET is mm/h£¬not mm/dt
Subroutine CREST_Excute(FDbasic, fDrain, fDpet, fDexport, RainStyle, PETstyle, &
		&DB, dT, Nmax, Nwarm, EchoLevel, conExport, ExpFormat, outN, expMean, expR, expC)
	Implicit None
	Character(*):: FDbasic, fDrain, fDpet, fDexport
	Integer(1):: RainStyle, PETstyle, EchoLevel, conExport(0:15), ExpFormat
	Integer(2):: DB(1:4)
	Integer(4):: Nmax, Nwarm, outN, expR(0:outN), expC(0:outN)
	Real(4):: dT
	Real(4):: expMean(1:(14+outN), 1:Nmax)
	
	Character(80):: fBasic(1:24), fileName, msg, str
	Integer(4):: I, nCols, nRows, lR, lC
	Real(4):: xllCor, yllCor, ceSize, noData
	Real(4), Allocatable:: temV(:,:)
	Logical:: fExist
	
	Real(4), Allocatable:: DEM(:,:), Dire(:,:), Cumu(:,:)
	Real(4), Allocatable:: WUin(:,:), WLin(:,:), WDin(:,:)
	Real(4), Allocatable:: StoSur(:,:), StoBas(:,:)
	Real(4), Allocatable:: WUmax(:,:), WLmax(:,:), WDmax(:,:)
	Real(4), Allocatable:: pIM(:,:), pB(:,:), pFc(:,:), pKe(:,:), pWcr(:,:), pEC(:,:)
	Real(4), Allocatable:: rSS(:,:), rSB(:,:), rSR(:,:), rTh(:,:), rSE(:,:), rGM(:,:), rLc(:,:)
	Real(4), Allocatable:: RPQ(:,:,:)

	Real(4):: lenSN, lenEW, LenCross, SpeedVegLocal, SpeedVegNext, TotalArea
	Integer(4):: TotalNum
	Logical:: SpecialPoint
	Real(4), Allocatable:: Area(:,:), Slope(:,:), AreaUpper(:,:)
	Real(4), Allocatable:: NextR(:,:), NextC(:,:), NextLen(:,:)
	Real(4), Allocatable:: SpeedSur(:,:), Speedbas(:,:)
	Real(4), Allocatable:: NextTimeSur(:,:), NextTimebas(:,:)
	
	Real(4), Allocatable:: toRsurA(:,:), toCsurA(:,:), toPsurA(:,:)
	Real(4), Allocatable:: toRsurB(:,:), toCsurB(:,:), toPsurB(:,:)
	Real(4), Allocatable:: toRbasA(:,:), toCbasA(:,:), toPbasA(:,:)
	Real(4), Allocatable:: toRbasB(:,:), toCbasB(:,:), toPbasB(:,:)

	Integer(4):: upL, N, L
	Integer(4), Allocatable:: posR(:), posC(:), outL(:)
	Real(4), Allocatable:: upA(:), StoOld(:), StoNew(:), ExcCur(:)

	Integer(4):: iT, toR, toC, Tb, Tx
	Integer(2):: Dcur(1:4), Tcur(1:3)
	Real(4):: toRoute
	
	Real(4), Allocatable:: WUout(:,:), WLout(:,:), WDout(:,:)
	Real(4), Allocatable:: Rain(:,:), PET(:,:), actE(:,:)
	Real(4), Allocatable:: ExcSur(:,:), ExcBas(:,:)
	Real(4), Allocatable:: RtdSur(:,:), RtdBas(:,:), DisQ(:,:)
	

	Call Nowmsg(msg)
	Open(199, file=Trim(fdExport) // "CREST." // Trim(msg) // ".log", form='formatted')! Log File
	If (EchoLevel>0) Then
		Write(*,"(A14, $)") "Basic Folder:"
		Write(*,*) Trim(fDbasic)
		Write(*,"(A14, $)") "Export Folder:"
		Write(*,*) Trim(fDexport)
		Write(*,"(A14, $)") "Rain Folder:"
		Write(*,*) Trim(fDrain)
		Write(*,"(A14, $)") "PET Folder:"
		Write(*,*) Trim(fDpet)
		Write(*,"(A14, I5, 3I3)") "Begin Data:", DB
		Write(*,"(A14, I5)") "Tatal Steps:", Nmax
		Write(*,"(A14, I5)") "Warm Up Steps:", Nwarm
		Write(*,"(A14, F5.2)") "dT:", dT
		Write(*,"(A14, I5)") "Control Point:", outN
		Do I=1, outN
			Write(*,"(A10, I3, A1, 2I5)") "Point", I, ":", expR(i), expC(i)
		End Do
		Write(*,"(A14, $)") "Routing Debug:"
		Write(*,"(I2)") conExport(0)
		Write(*,"(A14, $)") "Rain:"
		Write(*,"(I2)") conExport(1)
		Write(*,"(A14, $)") "PET:"
		Write(*,"(I2)") conExport(2)
		Write(*,"(A14, $)") "actE:"
		Write(*,"(I2)") conExport(3)
		Write(*,"(A14, $)") "ExcSur:"
		Write(*,"(I2)") conExport(4)
		Write(*,"(A14, $)") "ExcBas:"
		Write(*,"(I2)") conExport(5)
		Write(*,"(A14, $)") "Excess:"
		Write(*,"(I2)") conExport(6)
		Write(*,"(A14, $)") "StoSur:"
		Write(*,"(I2)") conExport(7)
		Write(*,"(A14, $)") "StoBas:"
		Write(*,"(I2)") conExport(8)
		Write(*,"(A14, $)") "Storage:"
		Write(*,"(I2)") conExport(9)
		Write(*,"(A14, $)") "WU:"
		Write(*,"(I2)") conExport(10)
		Write(*,"(A14, $)") "WL:"
		Write(*,"(I2)") conExport(11)
		Write(*,"(A14, $)") "WD:"
		Write(*,"(I2)") conExport(12)
		Write(*,"(A14, $)") "Discharge:"
		Write(*,"(I2)") conExport(13)
		Write(*,"(A14, $)") "For Next:"
		Write(*,"(I2)") conExport(14)
	Endif


	fBasic(1) = Trim(FDbasic) // "Dem.bif"
	fBasic(2) = Trim(FDbasic) // "Dire.bif"
	fBasic(3) = Trim(FDbasic) // "Cumu.bif"
	fBasic(4) = Trim(FDbasic) // "StoSur.bif"
	fBasic(5) = Trim(FDbasic) // "StoBas.bif"
	
	fBasic(6) = Trim(FDbasic) // "WUin.bif" ! Unit is %
	fBasic(7) = Trim(FDbasic) // "WLin.bif" ! Unit is %
	fBasic(8) = Trim(FDbasic) // "WDin.bif" ! Unit is %

	fBasic(9) = Trim(FDbasic) // "rSS.bif"
	fBasic(10) = Trim(FDbasic) // "rSB.bif"
	fBasic(11) = Trim(FDbasic) // "rSR.bif" ! Unit is %
	fBasic(12) = Trim(FDbasic) // "rTh.bif"
	fBasic(13) = Trim(FDbasic) // "rSE.bif"
	fBasic(14) = Trim(FDbasic) // "rGM.bif"
	fBasic(15) = Trim(FDbasic) // "rLc.bif"
	fBasic(16) = Trim(FDbasic) // "WUmax.bif"
	fBasic(17) = Trim(FDbasic) // "WLmax.bif"
	fBasic(18) = Trim(FDbasic) // "WDmax.bif"
	fBasic(19) = Trim(FDbasic) // "pIM.bif" ! Unit is %
	fBasic(20) = Trim(FDbasic) // "pB.bif"
	fBasic(21) = Trim(FDbasic) // "pFc.bif"
	fBasic(22) = Trim(FDbasic) // "pKe.bif"
	fBasic(23) = Trim(FDbasic) // "pWcr.bif"
	fBasic(24) = Trim(FDbasic) // "pEC.bif"
	Do I=1, 24
		Inquire(file=Trim(fBasic(I)), exist=fExist)
		If (.not. fExist) Then
			Write(*,"(A14, $)") "Missing File "
			Write(*,*) Trim(fBasic(I))
			Stop ! All files are necessary
		End If
	Enddo
	
	Call LoadBifHead(Trim(fBasic(1)), nCols, nRows, xLLcor, yLLcor, ceSize, noData)
	Allocate(DEM(1:nCols, 1:nRows), Dire(1:nCols, 1:nRows), Cumu(1:nCols, 1:nRows))
	Call LoadBif(fBasic(1), DEM, nCols, nRows)
	Call LoadBif(fBasic(2), Dire, nCols, nRows)
	Call LoadBif(fBasic(3), Cumu, nCols, nRows)
	
	Allocate(rSS(1:nCols, 1:nRows), rSB(1:nCols, 1:nRows), rSR(1:nCols, 1:nRows), rTh(1:nCols, 1:nRows))
	Call LoadBif(fBasic(9), rSS, nCols, nRows)
	Call LoadBif(fBasic(10), rSB, nCols, nRows)
	Call LoadBif(fBasic(11), rSR, nCols, nRows)
	Call LoadBif(fBasic(12), rTh, nCols, nRows)

	Allocate(rSE(1:nCols, 1:nRows), rGM(1:nCols, 1:nRows), rLc(1:nCols, 1:nRows))
	Call LoadBif(fBasic(13), rSE, nCols, nRows)
	Call LoadBif(fBasic(14), rGM, nCols, nRows)
	Call LoadBif(fBasic(15), rLc, nCols, nRows)

	Allocate(WUmax(1:nCols, 1:nRows), WLMax(1:nCols, 1:nRows), WDmax(1:nCols, 1:nRows))
	Call LoadBif(fBasic(16), WUmax, nCols, nRows)
	Call LoadBif(fBasic(17), WLmax, nCols, nRows)
	Call LoadBif(fBasic(18), WDmax, nCols, nRows)

	Allocate(pIM(1:nCols, 1:nRows), pB(1:nCols, 1:nRows), pFc(1:nCols, 1:nRows))
	Call LoadBif(fBasic(19), pIM, nCols, nRows)
	Call LoadBif(fBasic(20), pB, nCols, nRows)
	Call LoadBif(fBasic(21), pFc, nCols, nRows)

	Allocate(pKe(1:nCols, 1:nRows), pWcr(1:nCols, 1:nRows), pEC(1:nCols, 1:nRows))
	Call LoadBif(fBasic(22), pKe, nCols, nRows)
	Call LoadBif(fBasic(23), pWcr, nCols, nRows)
	Call LoadBif(fBasic(24), pEC, nCols, nRows)	
	
	! If found status file outputed from last time£¬read it as Initial value
	Call D2yyyymmddhh(DB, msg)
	Allocate(StoSur(1:nCols, 1:nRows), StoBas(1:nCols, 1:nRows))
	Allocate(WUin(1:nCols, 1:nRows), WLin(1:nCols, 1:nRows), WDin(1:nCols, 1:nRows))
	
	fileName=Trim(fDexport) // "CREST."// Trim(msg) // ".StoSur.bif"
	Inquire(File=Trim(fileName), Exist=fExist)
	If (.not. fExist) filename=fBasic(4)
	Call Loadbif(fileName, StoSur, nCols, nRows)

	fileName=Trim(fDexport) // "CREST."// Trim(msg) // ".StoBas.bif"
	Inquire(File=Trim(fileName), Exist=fExist)
	If (.not. fExist) filename=fBasic(5)
	Call Loadbif(fileName, StoBas, nCols, nRows)

	fileName=Trim(fDexport) // "CREST."// Trim(msg) // ".WU.bif"
	Inquire(File=Trim(fileName), Exist=fExist)
	If (.not. fExist) filename=fBasic(6)
	Call Loadbif(fileName, WUin, nCols, nRows)
	
	fileName=Trim(fDexport) // "CREST."// Trim(msg) // ".WL.bif"
	Inquire(File=Trim(fileName), Exist=fExist)
	If (.not. fExist) filename=fBasic(7)
	Call Loadbif(fileName, WLin, nCols, nRows)
	
	fileName=Trim(fDexport) // "CREST."// Trim(msg) // ".WD.bif"
	Inquire(File=Trim(fileName), Exist=fExist)
	If (.not. fExist) filename=fBasic(8)
	Call Loadbif(fileName, WDin, nCols, nRows)

	Do LR=1, nRows
		Do LC=1, nCols
			If (DEM(LC,LR)/=noData) Then
				WUin(LC,LR)=WUin(LC,LR)*WUmax(LC,LR)/100. ! Change unit of Soil Water from % to mm
				WLin(LC,LR)=WLin(LC,LR)*WLmax(LC,LR)/100.
				WDin(LC,LR)=WDin(LC,LR)*WDmax(LC,LR)/100.
				pFc(lC, lR) = pFc(lC, lR) * dT
			endif
		Enddo
	Enddo
	Write(199,"(A14, A4)") "Pretreat:", "0/4"
	If (EchoLevel>0) Write(*,"(A14, A4)") "Pretreat:", "0/4"

	!Load Grid based Threshold
	Allocate(temV(1:nCols, 1:nRows))	
	If (conExport(15) > 0) Then
		Allocate(RPQ(1:nCols, 1:nRows, 1:conExport(15)))	
		Do I = 1, conExport(15)
			Write(msg,"(I1)") I
			fileName = trim(FDbasic) // "RP" // trim(msg) // ".bif"
			Inquire(File=Trim(fileName), Exist=fExist)
            If (fExist) Then
				Call Loadbif(fileName, temV, nCols, nRows)
                Do lR = 1, nRows
                    Do lC = 1, nCols
                        RPQ(lC, lR, I) = temV(lC, lR)
                    End Do
                End Do
            Else
                DO lR = 1, nRows
                    Do lC = 1, nCols
                        RPQ(lC, lR, I) = noData
                    End Do
                End Do
            End If
        End Do
	End If


	! Routing pretreat
	Allocate(Area(1:nCols, 1:nRows), Slope(1:nCols, 1:nRows))
	Allocate(NextR(0:nCols, 0:nRows), NextC(0:nCols, 0:nRows))
	Allocate(NextLen(1:nCols, 1:nRows) )
	Allocate(SpeedSur(1:nCols, 1:nRows), Speedbas(1:nCols, 1:nRows))
	Allocate(NextTimeSur(0:nCols, 0:nRows), NextTimebas(0:nCols, 0:nRows))
	
	lenSN = ceSize * 120000 ! unit is M
	TotalNum = 0
	TotalArea = 0
	Do lR = 1, nRows
		Do lC = 1, nCols
			If (DEM(lC, lR) /= noData) Then
				lenEW = yllCor + (nRows - lR + 0.5) * ceSize
				lenEW = lenSN * Cos(lenEW / 180 * 3.1415926)
				LenCross = Sqrt((lenEW ** 2 + lenSN ** 2))
				Area(lC, lR) = lenSN * lenEW
				TotalArea = TotalArea + Area(lC, lR)
				TotalNum = TotalNum + 1
				
				SpecialPoint = .False.
				Select Case (Int(Dire(lC, lR)))
				Case (1)
					NextR(lC, lR) = lR - 1;  NextC(lC, lR) = lC;		NextLen(lC, lR) = lenSN
				Case (2)
					NextR(lC, lR) = lR - 1;  NextC(lC, lR) = lC + 1;	NextLen(lC, lR) = LenCross
				Case (3)
					NextR(lC, lR) = lR;	  NextC(lC, lR) = lC + 1;	NextLen(lC, lR) = lenEW
				Case (4)
					NextR(lC, lR) = lR + 1;  NextC(lC, lR) = lC + 1;	NextLen(lC, lR) = LenCross
				Case (5)
					NextR(lC, lR) = lR + 1;  NextC(lC, lR) = lC;		NextLen(lC, lR) = lenSN
				Case (6)
					NextR(lC, lR) = lR + 1;  NextC(lC, lR) = lC - 1;	NextLen(lC, lR) = LenCross
				Case (7)
					NextR(lC, lR) = lR;	  NextC(lC, lR) = lC - 1;	NextLen(lC, lR) = lenEW
				Case (8)
					NextR(lC, lR) = lR - 1;  NextC(lC, lR) = lC - 1;	NextLen(lC, lR) = LenCross
				Case Default
					NextR(lC, lR) = 0;		NextC(lC, lR) = 0;		NextLen(lC, lR) = lenSN 
				End Select
						
				If (NextR(lC, lR) < 1 .Or. nRows < NextR(lC, lR) .Or. NextC(lC, lR) < 1 .Or. nCols < NextC(lC, lR)) Then
					SpecialPoint = .True.
					NextR(lC, lR) = 0
					NextC(lC, lR) = 0
					NextLen(lC, lR) = lenSN
				Else
					If (DEM(NextC(lC, lR), NextR(lC, lR)) == noData) Then
						SpecialPoint = .True.
						NextR(lC, lR) = 0
						NextC(lC, lR) = 0
						NextLen(lC, lR) = lenSN
					End If
				End If

				SpeedVegLocal = 0.5 
				If (SpecialPoint) Then
					SpeedVegNext = SpeedVegLocal
					Slope(lC, lR) = rGM(lC, lR) / NextLen(lC, lR)
				Else
					SpeedVegNext = 0.5
					If (DEM(lC, lR) > DEM(NextC(lC, lR), NextR(lC, lR))) Then
						Slope(lC, lR) = (DEM(lC, lR) - DEM(NextC(lC, lR), NextR(lC, lR))) / NextLen(lC, lR)
					Else
						Slope(lC, lR) = rGM(lC, lR) / NextLen(lC, lR)
					End If
				End If
		
				SpeedSur(lC, lR) = rSS(lC, lR) * (SpeedVegLocal + SpeedVegNext) / 2 * (Slope(lC, lR)) ** rSE(lC, lR)
				If (Cumu(lC, lR) > rTh(lC, lR))  SpeedSur(lC, lR) = SpeedSur(lC, lR) * rSR(lC, lR)
				NextTimeSur(lC, lR) = NextLen(lC, lR) / SpeedSur(lC, lR) / 3600 

				Speedbas(lC, lR) = rSB(lC, lR) / 100 * rSS(lC, lR) * (SpeedVegLocal + SpeedVegNext) / 2 * (Slope(lC, lR)) ** rSE(lC, lR)
				If (Cumu(lC, lR) > rTh(lC, lR))  Speedbas(lC, lR) = Speedbas(lC, lR) * rSR(lC, lR)
				NextTimebas(lC, lR) = NextLen(lC, lR) / Speedbas(lC, lR) / 3600 
			End If
		End Do
	End Do
	NextTimeSur(0, 0) = 9999
	NextTimebas(0, 0) = 9999
	Write(199,"(A14, A4)") "Pretreat:", "1/4"
	If (EchoLevel>0) Write(*,"(A14, A4)") "Pretreat:", "1/4"


	If (conExport(0) > 0) Then
		Do lR = 1, nRows
			Do lC = 1, nCols
				If (DEM(lC, lR) == noData) Then
					Slope(lC, lR) = noData
					Area(lC, lR) = noData
				End If
			End Do
		End Do
		Call SaveFile(fDexport // "CREST." // ".Area", Area, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
		Call SaveFile(fDexport // "CREST." // ".Slope", Slope, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
	End If
	DeAllocate(rSS, rSB, rSR, rSE, rGM, SpeedSur, Speedbas, NextLen, Slope)


	Allocate(toRsurA(1:nCols, 1:nRows), toCsurA(1:nCols, 1:nRows), toPsurA(1:nCols, 1:nRows))
	Allocate(toRsurB(1:nCols, 1:nRows), toCsurB(1:nCols, 1:nRows), toPsurB(1:nCols, 1:nRows))
	Do lR = 1, nRows
		Do lC = 1, nCols
			If (DEM(lC, lR) /= noData) Then
				toRsurB(lC, lR) = lR
				toCsurB(lC, lR) = lC
				toPsurB(lC, lR) = 0
				Do While (toPsurB(lC, lR) < dT)
					toRsurA(lC, lR) = toRsurB(lC, lR)
					toCsurA(lC, lR) = toCsurB(lC, lR)
					toPsurA(lC, lR) = toPsurB(lC, lR)

					toRsurB(lC, lR) = NextR(toCsurA(lC, lR), toRsurA(lC, lR))
					toCsurB(lC, lR) = NextC(toCsurA(lC, lR), toRsurA(lC, lR))
					toPsurB(lC, lR) = toPsurB(lC, lR) + NextTimeSur(toCsurA(lC, lR), toRsurA(lC, lR))
				End Do
				toPsurB(lC, lR) = (dT - toPsurA(lC, lR)) / (toPsurB(lC, lR) - toPsurA(lC, lR))
				toPsurA(lC, lR) = 1 - toPsurB(lC, lR)
			End If
		End Do
	End Do
	Deallocate(NextTimeSur)
	Write(199,"(A14, A4)") "Pretreat:", "2/4"
	If (EchoLevel>0) Write(*,"(A14, A4)") "Pretreat:", "2/4"


	Allocate(toRBasA(1:nCols, 1:nRows), toCBasA(1:nCols, 1:nRows), toPBasA(1:nCols, 1:nRows))
	Allocate(toRBasB(1:nCols, 1:nRows), toCBasB(1:nCols, 1:nRows), toPBasB(1:nCols, 1:nRows))
	Do lR = 1, nRows
		Do lC = 1, nCols
			If (DEM(lC, lR) /= noData) Then
				toRBasB(lC, lR) = lR
				toCBasB(lC, lR) = lC
				toPBasB(lC, lR) = 0
				Do While (toPBasB(lC, lR) < dT)
					toRBasA(lC, lR) = toRBasB(lC, lR)
					toCBasA(lC, lR) = toCBasB(lC, lR)
					toPBasA(lC, lR) = toPBasB(lC, lR)

					toRBasB(lC, lR) = NextR(toCBasA(lC, lR), toRBasA(lC, lR))
					toCBasB(lC, lR) = NextC(toCBasA(lC, lR), toRBasA(lC, lR))
					toPBasB(lC, lR) = toPBasB(lC, lR) + NextTimeBas(toCBasA(lC, lR), toRBasA(lC, lR))
				End Do
				toPBasB(lC, lR) = (dT - toPBasA(lC, lR)) / (toPBasB(lC, lR) - toPBasA(lC, lR))
				toPBasA(lC, lR) = 1 - toPBasB(lC, lR)
			End If
		End Do
	End Do
	Deallocate(NextTimeBas)
	Write(199,"(A14, A4)") "Pretreat:", "3/4"
	If (EchoLevel>0) Write(*,"(A14, A4)") "Pretreat:", "3/4"
	

	Allocate(posR(1:TotalNum), posC(1:TotalNum), outL(1:TotalNum))
	Allocate(StoOld(1:TotalNum), StoNew(1:TotalNum), ExcCur(1:TotalNum), upA(1:TotalNum))
	
	upL = 1
	Do lR = 1, nRows
		Do lC = 1, nCols
			If (DEM(lC, lR) /= noData .And. NextR(lC, lR) == 0 .And. NextC(lC, lR) == 0) Then
				posR(upL) = lR
				posC(upL) = lC
				outL(upL) = Cumu(lC, lR)
				upL = upL + 1
			End If
		End Do
	End Do
	
	N = 1
	expC(0) = posC(N)
	expR(0) = posR(N)
	Do L = 1, upL-1
		If (outL(L) > outL(N)) Then
			N = L
			expC(0) = posC(N)
			expR(0) = posR(N)
		End If
	End Do
	Do L = 1, upL-1
		outL(L) = 0
	End Do
	Deallocate(NextR, NextC)
	

	L = 1
	Do While (L < upL)
		lR = posR(L)
		lC = posC(L)
		Call TestPixelDire(Dire, LR-1, LC, 5, nCols, nRows, L, upL, posR, posC, outL, TotalNum)
		Call TestPixelDire(Dire, LR-1, LC+1, 6, nCols, nRows, L, upL, posR, posC, outL, TotalNum)
		Call TestPixelDire(Dire, LR, LC+1, 7, nCols, nRows, L, upL, posR, posC, outL, TotalNum)
		Call TestPixelDire(Dire, LR+1, LC+1, 8, nCols, nRows, L, upL, posR, posC, outL, TotalNum)
		Call TestPixelDire(Dire, LR+1, LC, 1, nCols, nRows, L, upL, posR, posC, outL, TotalNum)
		Call TestPixelDire(Dire, LR+1, LC-1, 2, nCols, nRows, L, upL, posR, posC, outL, TotalNum)
		Call TestPixelDire(Dire, LR, LC-1, 3, nCols, nRows, L, upL, posR, posC, outL, TotalNum)
		Call TestPixelDire(Dire, LR-1, LC-1, 4, nCols, nRows, L, upL, posR, posC, outL, TotalNum)
		L = L + 1
	End Do
	Do L = 1, TotalNum
		StoNew(L) = (StoBas(posC(L), posR(L)) + StoSur(posC(L), posR(L))) * Area(posC(L), posR(L)) / 1000 
		upA(L) = Area(posC(L), posR(L))
	End Do
	Do L = TotalNum, 1, -1
		If (outL(L) /= 0) Then
			StoNew(outL(L)) = StoNew(outL(L)) + StoNew(L)
			upA(outL(L)) = upA(outL(L)) + upA(L)
		End If
	End Do

	Allocate(AreaUpper(1:nCols, 1:nRows))
	AreaUpper = noData
	Do L = 1, TotalNum
		AreaUpper(posC(L), posR(L)) = upA(L)
	End Do
	If (conExport(0) > 0) Call SaveFile(fDexport // "CREST." // ".Area", AreaUpper, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
	Deallocate(upA)

	If (EchoLevel>0) Write(*,"(A14, A4)") "Pretreat:", "4/4"


	! Main
	Allocate(WUout(1:nCols, 1:nRows), WLout(1:nCols, 1:nRows), WDout(1:nCols, 1:nRows))
	Allocate(Rain(1:nCols, 1:nRows), PET(1:nCols, 1:nRows), actE(1:nCols, 1:nRows))
	Allocate(ExcSur(1:nCols, 1:nRows), ExcBas(1:nCols, 1:nRows), DisQ(1:nCols, 1:nRows))
	Allocate(RtdSur(1:nCols, 1:nRows), RtdBas(1:nCols, 1:nRows))

	Call NowTime(Tcur)
	Tb=Tcur(1)*3600+Tcur(2)*60+Tcur(3)
	ExpMean=0
	Dcur=DB
	Do iT = 1, Nmax
		PET=0 
		Select Case (PETstyle)
		Case (0)
			Call D2mm(Dcur, msg)
			fileName = Trim(fDpet) // Trim(msg) // ".bif"
			Inquire(File=Trim(fileName), Exist=fExist)
			If (fExist) Then
				Call LoadBif(fileName, PET, nCols, nRows)
			End if
		Case (1)
			Call D2mm(Dcur, msg)
			fileName = Trim(fDpet) // "PET025." // Trim(msg) // ".bif"
			PET=0
			Inquire(File=Trim(fileName), Exist=fExist)
			If (fExist) Then
				Call LoadBIFbyRegion(Trim(fileName), PET, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
			Else
				Write(199,*) "Missing PET File: " // Trim(msg)
				If (EchoLevel>0) Write(*,*) "Missing PET File: " // Trim(msg)
			End if
		End Select
		
		Rain=0
		Select Case (RainStyle)
		Case (0)
			Call D2yyyymmddhh(Dcur, msg)
			fileName = Trim(fDrain) // Trim(msg) // ".bif"
			Inquire(File=Trim(fileName), Exist=fExist)
			If (fExist) Then
				Call LoadBif(fileName, Rain, nCols, nRows)
			End if
		Case (1)
			Call D2yyyymmddhh(Dcur, msg)
			fileName = Trim(fDrain) // "3B42RT." // Trim(msg) // ".6A.bin"
			Inquire(File=Trim(fileName), Exist=fExist)
			If (.not. fExist) Then
				fileName = Trim(fDrain) // "3B42RT." // Trim(msg) // ".6.bin"
				Inquire(File=Trim(fileName), Exist=fExist)
				If (.Not. fExist) Then
					fileName = Trim(fDrain) // "3B42RT." // Trim(msg) // ".bin"
					Inquire(File=Trim(fileName), Exist=fExist)
					If (.Not. fExist) Then
						Call D2yyyy(Dcur, str)
						fileName = Trim(fDrain) // Trim(str) // "/3B42RT." // Trim(msg) // ".bin"
						Inquire(File=Trim(fileName), Exist=fExist)
						If (.Not. fExist) Then
							fileName = Trim(fDrain) // Trim(str) // "/3B42RT." // Trim(msg) // ".6.bin"
							Inquire(File=Trim(fileName), Exist=fExist)
							If (.Not. fExist) fileName = Trim(fDrain) // Trim(str) // "/3B42RT." // Trim(msg) // ".6A.bin"
						End if
					End if
				End if
			End If

			Inquire(File=Trim(fileName), Exist=fExist)
			If (fExist) Then
				Call LoadRTByRegion(Trim(fileName), Rain, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
			Else
				Write(199,*) "Missing Rain File: " // Trim(msg)
				If (EchoLevel>0) Write(*,*) "Missing Rain File: " // Trim(msg)
			End if
		Case (11)
			Call D2yyyymmdd(Dcur, msg)
            fileName = trim(fDrain) // "RT." // Trim(msg) // ".Day.bif"
			Inquire(File=Trim(fileName), Exist=fExist)
            If (fExist) Then
				Call LoadBIFbyRegion(fileName, Rain, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
            Else
                Write(199,*) "Missing Rain File: " // Trim(msg)
				If (EchoLevel>0) Write(*,*) "Missing Rain File: " // Trim(msg)
            End If
		Case (99)
			Do lR = 1, nRows
				Do lC = 1, nCols
					If (DEM(lC, lR) /= noData) Rain(lc,lr)=10
				End Do
			End Do
		End Select	


		Do lR = 1, nRows
			Do lC = 1, nCols
				If (DEM(lC, lR) /= noData) Then
					If (PET(lC, lR) <= 0) Then
						PET(lC, lR) = 0
					Else
						Select Case (PETstyle)
						Case (1)
							PET(lC, lR) = PET(lC, lR) * dT * pKe(LC,LR)
						End Select
					End If
					
					If (Rain(lC, lR) <= 0) Then
						Rain(lC, lR) = 0
					Else
						Select Case (RainStyle)		
						Case (1)
							Rain(lC, lR) = Rain(lC, lR) * dT
						Case (11)
							!Rain(lC, lR) = Rain(lC, lR)
						End Select
					End If
				End if
			End Do
		End Do

		
		RtdSur=0
		RtdBas=0 
		Do lR = 1, nRows
			Do lC = 1, nCols
				If (DEM(lC, lR) /= noData) Then
					Call STOREFULL(WUin(LC,LR), WLin(LC,LR), WDin(LC,LR), WUout(LC,LR), WLout(LC,LR), WDout(LC,LR), &
						&Rain(LC,LR), PET(LC,LR), ExcSur(LC,LR), ExcBas(LC,LR), actE(LC,LR), WUmax(LC,LR), WLMax(LC,LR), WDmax(LC,LR), &
						&pIM(LC,LR), pB(LC,LR), pFc(LC,LR), pWcr(LC,LR), pEc(LC,LR))

					StoSur(lC, lR) = StoSur(lC, lR) + ExcSur(lC, lR)
					toRoute = StoSur(lC, lR) * rLc(lC, lR)
					StoSur(lC, lR) = StoSur(lC, lR) - toRoute
					toR = toRsurA(lC, lR)
					toC = toCsurA(lC, lR)
					If (toR > 0 .And. toC > 0) RtdSur(toC, toR) = RtdSur(toC, toR) + toRoute * toPsurA(lC, lR) * Area(lC, lR) / Area(toC, toR)
					toR = toRsurB(lC, lR)
					toC = toCsurB(lC, lR)
					If (toR > 0 .And. toC > 0) RtdSur(toC, toR) = RtdSur(toC, toR) + toRoute * toPsurB(lC, lR) * Area(lC, lR) / Area(toC, toR)
					
					StoBas(lC, lR) = StoBas(lC, lR) + ExcBas(lC, lR)
					toRoute = StoBas(lC, lR) * rLc(lC, lR)
					StoBas(lC, lR) = StoBas(lC, lR) - toRoute
					toR = toRbasA(lC, lR)
					toC = toCbasA(lC, lR)
					If (toR > 0 .And. toC > 0) RtdBas(toC, toR) = RtdBas(toC, toR) + toRoute * toPbasA(lC, lR) * Area(lC, lR) / Area(toC, toR)
					toR = toRbasB(lC, lR)
					toC = toCbasB(lC, lR)
					If (toR > 0 .And. toC > 0) RtdBas(toC, toR) = RtdBas(toC, toR) + toRoute * toPbasB(lC, lR) * Area(lC, lR) / Area(toC, toR)
				End If
			End Do
		End Do

		Do lR = 1, nRows 
			Do lC = 1, nCols
				If (DEM(lC, lR) /= noData) Then
					StoSur(lC, lR) = StoSur(lC, lR) + RtdSur(lC, lR)
					StoBas(lC, lR) = StoBas(lC, lR) + RtdBas(lC, lR)
				End If
			End Do
		End Do


		Do L = 1, TotalNum
			StoOld(L) = StoNew(L)
			StoNew(L) = 0
			ExcCur(L) = 0
		End Do
		Do L = 1, TotalNum
			StoNew(L) = (StoBas(posC(L), posR(L)) + StoSur(posC(L), posR(L))) * Area(posC(L), posR(L)) / 1000 
			ExcCur(L) = (ExcBas(posC(L), posR(L)) + ExcSur(posC(L), posR(L))) * Area(posC(L), posR(L)) / 1000
		End Do
		Do L = TotalNum, 1, -1
			If (outL(L) /= 0) Then
				StoNew(outL(L)) = StoNew(outL(L)) + StoNew(L)
				ExcCur(outL(L)) = ExcCur(outL(L)) + ExcCur(L)
			End If
		End Do
		Do L = 1, TotalNum 
			DisQ(posC(L), posR(L)) = (StoOld(L) + ExcCur(L) - StoNew(L)) / dT / 3600 !ÃëÁ¢Ã×
		End Do


		Do lR = 1, nRows
			Do lC = 1, nCols
				If (DEM(lC, lR) /= noData) Then
					If (Rain(lC, lR) > 0) expMean(1, iT) = expMean(1, iT) + Rain(lC, lR)
					If (PET(lC, lR) > 0) expMean(2, iT) = expMean(2, iT) + PET(lC, lR)
					If (actE(lC, lR) > 0) expMean(3, iT) = expMean(3, iT) + actE(lC, lR)
					
					If (ExcSur(lC, lR) > 0) expMean(4, iT) = expMean(4, iT) + ExcSur(lC, lR)
					If (ExcBas(lC, lR) > 0) expMean(5, iT) = expMean(5, iT) + ExcBas(lC, lR)
					
					If (StoSur(lC, lR) > 0) expMean(7, iT) = expMean(7, iT) + StoSur(lC, lR)
					If (StoBas(lC, lR) > 0) expMean(8, iT) = expMean(8, iT) + StoBas(lC, lR)
					
					WUin(lC, lR) = (WUin(lC, lR) + WUout(lC, lR)) / WUmax(lC, lR) * 50.
					WLin(lC, lR) = (WLin(lC, lR) + WLout(lC, lR)) / WLmax(lC, lR) * 50.
					WDin(lC, lR) = (WDin(lC, lR) + WDout(lC, lR)) / WDmax(lC, lR) * 50.
					If (WUin(lC, lR) > 0) expMean(10, iT) = expMean(10, iT) + WUin(lC, lR)
					If (WLin(lC, lR) > 0) expMean(11, iT) = expMean(11, iT) + WLin(lC, lR)
					If (WDin(lC, lR) > 0) expMean(12, iT) = expMean(12, iT) + WDin(lC, lR)
				End If
			End Do
		End Do

		expMean(1, iT) = expMean(1, iT) / TotalNum / dT
		expMean(2, iT) = expMean(2, iT) / TotalNum / dT
		expMean(3, iT) = expMean(3, iT) / TotalNum / dT
		
		expMean(4, iT) = expMean(4, iT) / TotalNum / dT
		expMean(5, iT) = expMean(5, iT) / TotalNum / dT
		expMean(6, iT) = expMean(4, iT) + expMean(5, iT)
		
		expMean(7, iT) = expMean(7, iT) / TotalNum / dT
		expMean(8, iT) = expMean(8, iT) / TotalNum / dT
		expMean(9, iT) = expMean(7, iT) + expMean(8, iT)
		
		expMean(10, iT) = expMean(10, iT) / TotalNum
		expMean(11, iT) = expMean(11, iT) / TotalNum
		expMean(12, iT) = expMean(12, iT) / TotalNum
		
		If (iT > 1) expMean(13, iT) = expMean(9, iT - 1) + expMean(6, iT) - expMean(9, iT)
		Do I = 0, outN
			expMean(14 + I, iT) = DisQ(expC(I), expR(I))
		End Do


		If (iT > Nwarm) Then
			Do lR = 1, nRows
				Do lC = 1, nCols
					If (DEM(lC, lR) == noData) Then
						Rain(lC, lR) = noData
						PET(lC, lR) = noData
						actE(lC, lR) = noData
						ExcSur(lC, lR) = noData
						ExcBas(lC, lR) = noData
						temV(lC, lR) = noData
						StoSur(lC, lR) = noData
						StoBas(lC, lR) = noData
						DisQ(lC, lR) = noData
					End If
				End Do
			End Do
			If (conExport(1) > 0) Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) // ".Rain", Rain, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			If (conExport(2) > 0) Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) // ".PET", PET, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			If (conExport(3) > 0) Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) // ".actE", actE, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			
			If (conExport(4) > 0) Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) // ".ExcSur", ExcSur, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			If (conExport(5) > 0) Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) // ".ExcBas", ExcBas, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			If (conExport(6) > 0) Then
				Do lR = 1, nRows
					Do lC = 1, nCols
						If (DEM(lC, lR) /= noData) temV(lC, lR) = ExcSur(lC, lR) + ExcBas(lC, lR)
					End Do
				End Do
				Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) // ".Excess", temV, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			End If
			
			If (conExport(7) > 0) Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) //  ".StoSur", StoSur, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			If (conExport(8) > 0) Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) //  ".StoBas", StoBas, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			If (conExport(9) > 0) Then
				Do lR = 1, nRows
					Do lC = 1, nCols
						If (DEM(lC, lR) /= noData) temV(lC, lR) = StoSur(lC, lR) + StoBas(lC, lR)
					End Do
				End Do
				Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) // ".Storage", temV, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			End If
			
			If (conExport(10) > 0) Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) //  ".WU", WUin, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			If (conExport(11) > 0) Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) //  ".WL", WLin, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			If (conExport(12) > 0) Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) //  ".WD", WDin, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			
			If (conExport(13) > 0) Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) //  ".Q", DisQ, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)

			If (conExport(15) > 0) Then
				Do lR = 1, nRows
                    Do lC = 1, nCols
                        If (DEM(lC, lR) == noData) Then
                            temV(lC, lR) = noData
                        Else
                            If (Cumu(lC, lR) >= rTh(lC, lR) * 10) Then
                                temV(lC, lR) = -1
                            Else
                                temV(lC, lR) = 0
                            End If
                            
                            Do I = 1, conExport(15)
                                If (DisQ(lC, lR) > RPQ(lC, lR, I)) temV(lC, lR) = I
                            End Do
                        End If
                    End Do
                End Do
				Call SaveFile(Trim(fDexport) // "CREST." // Trim(msg) //  ".Qlevel", temV, nCols, nRows, xLLcor, yLLcor, ceSize, noData, ExpFormat)
			End If
		End If
		
		
		Write(199,"(A14, I5, 3I3, $)") "Calculate:", Dcur
		If (EchoLevel>0) Then
			Write(*,"(A14, I5, 3I3, $)") "Calculate:", Dcur
			Call NowTime(Tcur)
			Tx = Tcur(1) * 3600 + Tcur(2) * 60 + Tcur(3)
			Tx = Real(Tx - Tb) / iT * (Nmax - iT)
			Tcur(1) = Tx / 3600
			Tx = Tx - Tcur(1) * 3600
			Tcur(2) = Tx / 60
			Tx = Tx - Tcur(2) * 60
			Tcur(3) = Tx
			Write(*,"(A20, 3I3)") "Time Left:", Tcur
		End if
		Write(199,"(A20, 3I3)") "Time Left:", Tcur


		WUin = WUout
		WLin = WLout
		WDin = WDout
		Call Date_Change(Dcur, 4, dT)
	End Do
	

	If (conExport(14) > 0) Then
		Call D2yyyymmddhh(Dcur, msg)
		Call SaveBif(Trim(fDexport) // "CREST." // Trim(msg) // ".StoSur.bif", WDout, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
		Call SaveBif(Trim(fDexport) // "CREST." // Trim(msg) // ".StoBas.bif", WDout, nCols, nRows, xLLcor, yLLcor, ceSize, noData)

		Do lR = 1, nRows
			Do lC = 1, nCols
				If (DEM(lC, lR) /= noData) Then
					WUout(lC, lR) = WUout(lC, lR) / WUmax(lC, lR) * 100.
					WLout(lC, lR) = WLout(lC, lR) / WLmax(lC, lR) * 100.
					WDout(lC, lR) = WDout(lC, lR) / WDmax(lC, lR) * 100.
				End If
			End Do
		End Do 	  
		Call SaveBif(Trim(fDexport) // "CREST." // Trim(msg) // ".WU.bif", WUout, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
		Call SaveBif(Trim(fDexport) // "CREST." // Trim(msg) // ".WL.bif", WLout, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
		Call SaveBif(Trim(fDexport) // "CREST." // Trim(msg) // ".WD.bif", WDout, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
	End If
	

		Write(199,"(A14, $)") "Basic Folder:"
		Write(199,*) Trim(fDbasic)
		Write(199,"(A14, $)") "Export Folder:"
		Write(199,*) Trim(fDexport)
		Write(199,"(A14, $)") "Rain Folder:"
		Write(199,*) Trim(fDrain)
		Write(199,"(A14, $)") "PET Folder:"
		Write(199,*) Trim(fDpet)
		Write(199,"(A14, I5, 3I3)") "Begin Data:", DB
		Write(199,"(A14, I5)") "Tatal Steps:", Nmax
		Write(199,"(A14, I5)") "Warm Up Steps:", Nwarm
		Write(199,"(A14, F5.2)") "dT:", dT
		Write(199,"(A14, I5)") "Control Point:", outN
		Do I=1, outN
			Write(199,"(A10, I3, A1, 2I5, F10.0, F15.3)") "Point", I, ":", expR(i), expC(i), Cumu(expC(I), expR(I)), AreaUpper(expC(I), expR(I))/1000000
		End Do
		Write(199,"(A14, $)") "Routing Debug:"
		Write(199,"(I2)") conExport(0)
		Write(199,"(A14, $)") "Rain:"
		Write(199,"(I2)") conExport(1)
		Write(199,"(A14, $)") "PET:"
		Write(199,"(I2)") conExport(2)
		Write(199,"(A14, $)") "actE:"
		Write(199,"(I2)") conExport(3)
		Write(199,"(A14, $)") "ExcSur:"
		Write(199,"(I2)") conExport(4)
		Write(199,"(A14, $)") "ExcBas:"
		Write(199,"(I2)") conExport(5)
		Write(199,"(A14, $)") "Excess:"
		Write(199,"(I2)") conExport(6)
		Write(199,"(A14, $)") "StoSur:"
		Write(199,"(I2)") conExport(7)
		Write(199,"(A14, $)") "StoBas:"
		Write(199,"(I2)") conExport(8)
		Write(199,"(A14, $)") "Storage:"
		Write(199,"(I2)") conExport(9)
		Write(199,"(A14, $)") "WU:"
		Write(199,"(I2)") conExport(10)
		Write(199,"(A14, $)") "WL:"
		Write(199,"(I2)") conExport(11)
		Write(199,"(A14, $)") "WD:"
		Write(199,"(I2)") conExport(12)
		Write(199,"(A14, $)") "Discharge:"
		Write(199,"(I2)") conExport(13)
		Write(199,"(A14, $)") "For Next:"
		Write(199,"(I2)") conExport(14)
	Close (199)

	If (EchoLevel>0) Then
		Write(*,"(A14, $)") "Finished:"
		Write(*,*) "CREST, Jiahu Wang, Jul.1 09"
	End if
End Subroutine




!===========================================================!
!	Written by Jiahu Wang
!	Excess Stirage Model, three layers, two water source
!===========================================================!

Subroutine STOREFULL(WUin, WLin, WDin, WUout, WLout, WDout, &
			&P, E, RS, Rb, EA, WUmax, WLMax, WDmax, &
			&pIM, pB, pFc, pWcr, pEc)
	Implicit None
	Real*4 WUin, WLin, WDin, WUout, WLout, WDout
	Real*4 P, E, Rs, Rb, EA, WUmax, WLMax, WDmax
	Real*4 pIM, pB, pFc, pKe, pWcr, pEc
	Real*4 Wm, Wmm, W, a, R, dW, pBtem

	pEc=0.0; pEc=0.2
	RS = 0
	Rb = 0
	If (P - E > 0) Then
		Wm = WUmax + WLMax + WDmax
		W = WUin + WLin + WDin
		pBtem = 1 + pB
		Wmm = pBtem * Wm
		a = Wmm * (1 - (1 - W / Wm) ** (1 / pBtem))
		If (P - E + a < Wmm) Then
			r = P - E - Wm * ((1 - a / Wmm) ** pBtem - (1 - (a + P - E) / Wmm) ** pBtem)
			If (r<0) r=0 
		Else
			r = P - E - (Wm - W)
		End If

		dW = P - E - r
		If (WUin + dW <= WUmax) Then
			WUout = WUin + dW
			WLout = WLin
			WDout = WDin
		Else
			WUout = WUmax
			dW = dW - (WUout - WUin)
			If (WLin + dW <= WLMax) Then
				WLout = WLin + dW
				WDout = WDin
			Else
				WLout = WLMax
				dW = dW - (WLout - WLin)
				If (WDin + dW <= WDmax) Then
					WDout = WDin + dW
				Else
					WDout = WDmax
					r = r + (WDin + dW - WDmax)
				End If
			End If
		End If
		
		r = r * (1 - pIM)
		If (P < pFc) Then
			Rb = r
			RS = 0
		Else
			Rb = r / P * pFc
			RS = r - Rb
		End If
		RS = RS + (P - E) * pIM
	Else
		dW = E - P
		If (WUin - dW > pWcr) Then
			WUout = WUin - dW
			WLout = WLin
			WDout = WDin
		Else
			WUout = pWcr
			dW = (dW - (WUin - WUout)) * pEc
			If (WLin - dW > 0) Then
				WLout = WLin - dW
				WDout = WDin
			Else
				WLout = pWcr
				dW = dW - (WLin - WLout)
				dW = dW * (WDin - pWcr) / (WDmax - pWcr)
				If (WDin - dW > pWcr) Then
					WDout = WDin - dW
				Else
					WDout = pWcr
				End If
			End If
		End If
	End If
	EA=Wuin+WLin+WDin-WUout-WLout-WDout+P-RS-Rb
End Subroutine


Subroutine LoadRTByRegion(filePath, Vout, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
	Implicit None
	Character(*):: filePath
	Integer*4:: nCols, nRows, souCols, souRows
	Real*4:: xLLcor, yLLcor, ceSize, noData, Vout(1:nCols, 1:nRows), souXll, souYll, souCS, souND
	Real*4, Allocatable:: Vsou(:,:)

	souRows = 480
	souCols = 1440
	souXll = -180
	souYll = -60
	souND = -9999.
	souCS = 0.25
	Allocate(Vsou(1:souCols, 1:souRows))
	Call LoadRT(filePath, Vsou, souCols, souRows, souND)
	Call GetBfromA(souCols, souRows, souXll, souYll, souCS, souND, Vsou, &
			&nCols, nRows, xllCor, yllCor, ceSize, noData, Vout)
End Subroutine


Subroutine LoadRT(filePath, Vsou, nCols, nRows, noData)
	Implicit None
	Character(*):: filePath
	Character*2880:: EmptySpace
	Integer*4:: nCols, nRows, LR, LC, BlockLength, gC
	Real*4:: noData, Vsou(1:nCols, 1:nRows)
	
	Integer*2, Allocatable:: Idir(:,:)
	Real*4, Allocatable:: temS(:,:)
	Allocate(Idir(1:nCols, 1:nRows), temS(1:nCols, 1:nRows))

	BlockLength=2880+nCols*nRows*2
	Open(102, file=Trim(filePath), access='direct', form='unformatted', recl=BlockLength)
		Read(102, rec=1)EmptySpace, ((Idir(LC,LR), LC=1, nCols), LR=1, nRows)
	Close(102)

	Do LR=1,nRows
		Do LC=1,nCols
			Call SwapInt(Idir(LC,LR))
		enddo
	enddo

	temS=noData
	Do LR = 1, nRows
		Do LC = 1, nCols
			If(Idir(LC, LR)>=0)Then
				temS(LC, LR) = Real(Idir(LC, LR)) / 100.0
			elseif (Idir(LC, LR)==-1)Then
				temS(LC, LR) = 0.0
			End If
		Enddo
	Enddo
	
	Do LC = 1, nCols
		If (LC <= 720) Then
			gC = LC + 720
		Else
			gC = LC - 720
		End If
		Do LR = 1, nRows
			Vsou(gC, LR) = temS(LC, LR)
		enddo
	enddo
End 


Subroutine SwapInt(Vsou)
	Implicit None
	Integer*2:: Vsou, temI
	Character*1:: cha2(1:2), cha
	Equivalence (temI, cha2)
	temI=Vsou
	cha=cha2(1)
	cha2(1)=cha2(2)
	cha2(2)=cha
	Vsou=temI
End 


Subroutine GetBfromA(AnCols, AnRows, AxLLcor, AyLLcor, AcellSize, AnoData, mtxA, &
					&BnCols, BnRows, BxLLcor, ByLLcor, BcellSize, BnoData, mtxB)
	Implicit None
	Integer*4:: AnCols, AnRows, BnCols, BnRows, LR, LC, souR, souC
	Real*4:: AxLLcor, AyLLcor, AcellSize, AnoData, mtxA(1:AnCols, 1:AnRows), &
			&BxLLcor, ByLLcor, BcellSize, BnoData, mtxB(1:BnCols, 1:BnRows), SR, SC

	mtxB=BnoData
	Do LR=1, BnRows
		SR = ByLLcor + BnRows * BcellSize - (LR - 0.5) * BcellSize
		souR = (AyLLcor + AnRows * AcellSize - SR) / AcellSize + 1
		If (souR>=1 .And. souR <= AnRows) Then
			Do LC=1, BnCols
				sC = (LC - 0.5) * BcellSize + BxLLcor
				souC = (sC - AxLLcor) / AcellSize + 1
				If (souC>=1 .And. souC <= AnCols) Then
					If (mtxA(souC, souR) /= AnoData)  mtxB(LC, LR) = mtxA(souC, souR)
				End If
			Enddo
		End If
	Enddo
End Subroutine


Subroutine Date_Change(D, Mark, Num)
	Integer*2:: D(1:4), Mark, I, Days
	Real(4):: Num
	
	D(Mark)=D(Mark)+Num
	If (d(4)>=24) Then ! Ê±
		i=d(4)/24
		d(3)=d(3)+i
		d(4)=d(4)-i*24
	End If
	
	If (d(2)==2) Then
		If (Mod(d(1), 4)==0) Then
			Days=29
		Else
			Days=28
		End If
	elseif (d(2)==4 .or. d(2)==6 .or. d(2)==9 .or. d(2)==11) Then
		Days=30
	Else
		Days=31
	End If
	If (d(3)>Days) Then
		d(3)=d(3)-Days
		d(2)=d(2)+1
	End If
	If (d(2)>12) Then
		d(2)=d(2)-12
		d(1)=d(1)+1
	End If
End


Subroutine D_yyyymmddhh(D, msg)
	Character(*):: msg
	Integer*2:: D(1:4)
	Character(4):: temA
	temA=msg(1:4);	 Read(temA,*) d(1)
	temA=msg(5:6);	 Read(temA,*) d(2)
	temA=msg(7:8);	 Read(temA,*) d(3)
	temA=msg(9:10);	 Read(temA,*) d(4)
End


Subroutine D2yyyymmdd(D, msg)
	Character(*):: msg
	Integer*2:: D(1:4), PETStyle
	Character(4):: temA

	msg=""
	Write(temA, "(I4)") d(1);									msg(1:4)=temA
	Write(temA, "(I2)") d(2);	If (d(2)<10) temA(1:1)="0";		msg(5:6)=temA
	Write(temA, "(I2)") d(3);	If (d(3)<10) temA(1:1)="0";		msg(7:8)=temA
End


Subroutine D2yyyymmddhh(D, msg)
	Character(*):: msg
	Integer*2:: D(1:4), PETStyle
	Character(4):: temA

	msg=""
	Write(temA, "(I4)") d(1);									msg(1:4)=temA
	Write(temA, "(I2)") d(2);	If (d(2)<10) temA(1:1)="0";		msg(5:6)=temA
	Write(temA, "(I2)") d(3);	If (d(3)<10) temA(1:1)="0";		msg(7:8)=temA
	Write(temA, "(I2)") d(4);	If (d(4)<10) temA(1:1)="0";		msg(9:10)=temA
End


Subroutine D2mm(D, msg)
	Character(*):: msg
	Integer*2:: D(1:4)

	msg=""
	Write(msg, "(I2)") d(2)
	If (d(2)<10) msg(1:1)="0"
End


Subroutine D2yyyy(D, msg)
	Character(*):: msg
	Integer*2:: D(1:4)

	msg=""
	Write(msg, "(I4)") d(1)
End


Subroutine Nowmsg(msg)
	Implicit None
	Character(*) msg
	Character(20) A, B
	Call Date_and_Time(A, B)
	msg=""
	msg(1:8)=A(1:8)
	msg(9:14)=B(1:6)
End Subroutine


Subroutine NowTime(T)
	Implicit None
	Character(80):: msg, temA
	Integer(2):: T(1:3)
	Call Time(msg)
	temA=msg(1:2);	 Read(temA,*) T(1)
	temA=msg(4:5);	 Read(temA,*) T(2)
	temA=msg(7:8);	 Read(temA,*) T(3)
End Subroutine


Subroutine RT_Date2Name(fileName, D)
	Character(*):: fileName
	Integer*2:: D(1:4)
	Character(4):: temA

	filename="3B42RT."
	Write(temA, "(I4)") d(1);	filename(8:11)=temA
	Write(temA, "(I2)") d(2);	If (d(2)<10) temA(1:1)="0";		filename(12:13)=temA
	Write(temA, "(I2)") d(3);	If (d(3)<10) temA(1:1)="0";		filename(14:15)=temA
	Write(temA, "(I2)") d(4);	If (d(4)<10) temA(1:1)="0";		filename(16:17)=temA
End


Subroutine LoadbifHead(filePath, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
	Implicit None
	Character(*):: filePath
	Integer*4:: nCols, nRows, LR, LC, BlockLength
	Real*4:: xLLcor, yLLcor, ceSize, noData
	  
	BlockLength=24
	Open (102, file=Trim(filePath), access='direct', form='unformatted', recl=BlockLength)
		Read(102, rec=1) nCols, nRows, xLLcor, yLLcor, ceSize, noData
	Close (102)
End Subroutine


Subroutine Loadbif(filePath, Vout, nCols, nRows)
	Implicit None
	Character(*):: filePath
	Character*26:: EmptySpace
	Integer*4:: nCols, nRows, LR, LC, BlockLength
	Real*4:: xLLcor, yLLcor, ceSize, noData, Vout(1:nCols,1:nRows)

	BlockLength=50+nRows* nCols
	Open (102, file=Trim(filePath), access='direct', form='unformatted', recl=BlockLength)
		Read(102, rec=1) nCols, nRows, xLLcor, yLLcor, ceSize, noData, EmptySpace, ((Vout(LC,LR), LC=1, nCols), LR=1, nRows) 
	Close (102)
End Subroutine


Subroutine Savebif(filePath, Vout, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
	Implicit None
	Character(*):: filePath
	Character*26:: EmptySpace
	Integer*4:: nCols, nRows, LR, LC, BlockLength
	Real*4:: xLLcor, yLLcor, ceSize, noData, Vout(1:nCols, 1:nRows)

	BlockLength=50+nRows*nCols
	Open(102, file=Trim(filePath), status="REPLACE", access='direct', form='unformatted', recl=BlockLength)
		Write(102, rec=1)nCols, nRows, xLLcor, yLLcor, ceSize, noData, EmptySpace, ((Vout(LC,LR), LC=1, nCols), LR=1, nRows) 
	Close (102)
End Subroutine


Subroutine LoadBIFbyRegion(filePath, Vout, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
	Implicit None
	Character(*):: filePath
	Integer*4:: nCols, nRows, souCols, souRows
	Real*4:: xLLcor, yLLcor, ceSize, noData, Vout(1:nCols, 1:nRows), souXll, souYll, souCS, souND
	Real*4, Allocatable:: Vsou(:,:)
	
	Call LoadbifHead(filePath, souCols, souRows, souXLL, souYLL, souCS, souND)
	Allocate(Vsou(1:souCols, 1:souRows))
	Call Loadbif(filePath, Vsou, souCols, souRows)
	Call GetBfromA(souCols, souRows, souXll, souYll, souCS, souND, Vsou, &
			&nCols, nRows, xllCor, yllCor, ceSize, noData, Vout)
End Subroutine


Subroutine TestPixelDire(Dire, LR, LC, GoalDire, nCols, nRows, L, upL, posR, posC, outL, TotalNum)
	Real(4):: Dire(1:nCols, 1:nRows)
	Integer(4):: LR, LC, nCols, nRows, GoalDire, L, upL, TotalNum, posR(1:TotalNum), posC(1:TotalNum), outL(1:TotalNum)
	If (1 <= LR .And. LR <= nRows .And. 1 <= LC .And. LC <= nCols) Then
		If (Int(Dire(LC, LR)) == GoalDire) Then
			posR(upL) = LR
			posC(upL) = LC
			outL(upL) = L
			upL=upL+1
		End If
	End If
End


Subroutine SaveFile(fileName, Vout, nCols, nRows, xLLcor, yLLcor, ceSize, noData, FileStyle)
	Implicit None
	Character(*):: fileName
	Integer(1):: FileStyle
	Integer(4):: nCols, nRows
	Real(4):: Vout(1:nCols, 1:nRows), xLLcor, yLLcor, ceSize, noData
	
	If (FileStyle==0) Then
		Call Savebif(Trim(fileName) // ".bif", Vout, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
	Else
		Call SaveAscii(Trim(fileName) // ".txt", Vout, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
	Endif
End Subroutine


Subroutine SaveAscii(fileName, Vout, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
	Implicit None
	Character(*):: fileName
	Integer(4):: lC, lR, nCols, nRows
	Real(4):: Vout(1:nCols, 1:nRows), xLLcor, yLLcor, ceSize, noData
	
	Open(102,file=Trim(fileName),form='formatted')
		Write(102,"('ncols 		')",advance='no'); Write(102,"(i8)")nCols
		Write(102,"('nrows 		')",advance='no'); Write(102,"(i8)")nRows
		Write(102,"('xllcorner 	')",advance='no'); Write(102,"(f10.6)")xllcor
		Write(102,"('yllcorner 	')",advance='no'); Write(102,"(f10.6)")yLLcor
		Write(102,"('cellsize 	 ')",advance='no'); Write(102,"(f10.6)")ceSize
		Write(102,"('NODATA_value ')",advance='no'); Write(102,"(f10.0)")noData
		Do lR=1, nRows
			Do lC=1, nCols
				Write(102,'(f10.3,A)',advance='no')Vout(lC,lR),' '
			enddo
			Write(102,'')
		enddo
	Close (102)
End Subroutine


Subroutine LoadASciiHead(fileName, nCols, nRows, xLLcor, yLLcor, ceSize, noData)
	Implicit None
	Character(*):: fileName
	Character(80):: msg
	Integer(4):: nCols, nRows
	Real(4):: xLLcor, yLLcor, ceSize, noData
	
	Open(102,file=Trim(fileName),form='formatted')
		Read(102,*) msg, nCols
		Read(102,*) msg, nRows
		Read(102,*) msg, xllcor
		Read(102,*) msg, yLLcor
		Read(102,*) msg, ceSize
		Read(102,*) msg, noData
	Close (102)
End Subroutine


Subroutine LoadAscii(fileName, Vout, nCols, nRows)
	Implicit None
	Character(*):: fileName
	Character(11):: msg
	Integer(4):: lC, lR, nCols, nRows
	Real(4):: Vout(1:nCols, 1:nRows), X
	
	Open(102,file=Trim(fileName),form='formatted')
		Read(102,*) msg, X
		Read(102,*) msg, X
		Read(102,*) msg, X
		Read(102,*) msg, X
		Read(102,*) msg, X
		Read(102,*) msg, X
		Do lR=1, nRows
			Do lC=1, nCols
				Read(102,*)Vout(lC,lR)
			enddo
		enddo
	Close (102)
End Subroutine