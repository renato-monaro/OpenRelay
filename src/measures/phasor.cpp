/*
# Copyright (C) 2008-2015 Rodolfo Varraschim Rocha
# This file is part of OpenRelay.
#
# OpenRelay is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenRelay is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
#
# RDFT is based on: A. Phadke, J. Thorp, M. Adamiak, A new measurement technique for
#   tracking phasor, local system frequency, and rate of change of frequency,
#   Power Apparatus and Systems, IEEE Transactions on 95 (5) (1983)
#   1025{1038. doi:10.1109/T-PAS.1976.32105.
#
# LS is based on: M. Sachdev, M. Baribeau, A new algorithm for digital impedance relays,
#   Power Apparatus and Systems, IEEE Transactions on PAS-98 (6) (1979)
#   2232{2240. doi:10.1109/TPAS.1979.319422.
#
# PLL is based on: A. Ziarani, A. Konrad, A novel method of extraction of nonstationary
#   sinusoids, in: Signal Processing 84 ELSEVIER, Vol. 1, 2004, pp. 1323{
#   1346. doi:10.1016/j.sigpro.2004.05.008.
#
# RLS is based on: I. Sadinezhad, V. Agelidis, Real-time power system phasors and harmon-
#   ics estimation using a new decoupled recursive-least-squares technique
#   for dsp implementation, Industrial Electronics, IEEE Transactions on
#   60 (6) (2013) 2295{2308. doi:10.1109/TIE.2012.2192895.
#
# DWT is based on: S. Brahma, R. Kavasseri, Investigating the performance of discrete
#   wavelet transform for phasor estimation in digital relays, in: Power
#   and Energy Society General Meeting - Conversion and Delivery of
#   Electrical Energy in the 21st Century, 2008 IEEE, 2008, pp. 1{6.
#   doi:10.1109/PES.2008.4596113.
#
# RWT is based on: J. Ren, M. Kezunovic, Elimination of dc offset in accurate phasor es-
#   timation using recursive wavelet transform, in: PowerTech, 2009 IEEE
#   Bucharest, 2009, pp. 1{5. doi:10.1109/PTC.2009.5282067.
#
*/

#include "phasor.h"

using namespace orelay;

RDFT::RDFT(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq, unsigned Harm){
        FunctionalDescription_set("Recursive Discrete Fourier Transform");
		Input=CH_IN;
		Output=CH_OUT;
		Harmonic=Harm;
		SystemFrequency=Freq;
		j=0;
		}
RDFT::RDFT(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq){
        Input=CH_IN;
		Output=CH_OUT;
		Harmonic=1;
		SystemFrequency=Freq;
		j=0;
		}
bool RDFT::Prepare(float Samp){
        RelaySamplingRate=Samp;
		N=(unsigned)(RelaySamplingRate/SystemFrequency);
		cout<<"Recursive Fourier N"<<N<<endl;
		return true;
		}
void RDFT::Run(){
        Complex F,i(0.0,1.0);
		F= Complex(0.0,0.0);

		if(j==0){
			for(unsigned p=0;p<N;p++)
			F+=Input->get_Value(p)*exp(((analogic)(2*M_PI/N)*i)*(analogic)(p*Harmonic));
			F/=(analogic)N;
			F*=sqrt(2);
			}
		else {
			F = Output->get_Value() + i*(2/(sqrt(2)*N))*(Input->get_Value()-Input->get_Value(N))*exp(i*(analogic)(2*M_PI/N)*(analogic)(Harmonic));
			}

		Output->insert_Value(F);

		//cout<<abs(Output->get_Value())<<" "<<arg(Output->get_Value())<<" "<<real(Output->get_Value())<<" "<<imag(Output->get_Value())<<endl;
		j++;
		if(j>=N)
			j=0;
		}

LS::LS(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq, unsigned Harm, unsigned HarmRange){
		    FunctionalDescription_set("Least Squares");
            Input=CH_IN;
            Output=CH_OUT;
            Harmonic=Harm;
            SystemFrequency=Freq;
            SpecRange=HarmRange;
            k=0;
		}
LS::LS(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq, unsigned Harm){
		    Input=CH_IN;
            Output=CH_OUT;
            Harmonic=Harm;
            SystemFrequency=Freq;
            SpecRange=5;
            k=0;
		}
LS::LS(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq){
		    Input=CH_IN;
            Output=CH_OUT;
            Harmonic=1;
            SystemFrequency=Freq;
            SpecRange=2;
            k=0;
		}
bool LS::Prepare(float Samp){
            RelaySamplingRate=Samp;
            N=(unsigned)(RelaySamplingRate/SystemFrequency);
            M= SpecRange;    //Numero de Harmônicos a se analisar
            cout<<"Least Squares N"<<N<<endl;
            return true;
        }
void LS::Run(){
            Complex Sum,i(0.0,1.0);
            Sum= Complex(0.0,0.0);

            gsl_matrix *vX = gsl_matrix_calloc(2*M+2, 1);
            gsl_matrix *jan = gsl_matrix_calloc(1, N);
            gsl_matrix *A = gsl_matrix_calloc(N, 2*M+2);
            gsl_matrix *AtA = gsl_matrix_calloc(2*M+2, 2*M+2);
            gsl_matrix *tmp = gsl_matrix_calloc(2*M+2, 2*M+2);
            gsl_matrix *inv_AtA = gsl_matrix_calloc(2*M+2, 2*M+2);
            gsl_matrix *PA = gsl_matrix_calloc(2*M+2, N);
            gsl_permutation *p = gsl_permutation_alloc(2*M+2);

            for(unsigned j=0;j<N;j++){
                nsamp=Input->get_Value((N-1)-j);
                gsl_matrix_set (jan, 0, j, nsamp);
                gsl_matrix_set (A, j, 0, 1);
                gsl_matrix_set (A, j, 1, (j+k)*(1/RelaySamplingRate));
                col=2;
                for (unsigned p=0;p<M;p++){
                    gsl_matrix_set (A, j, col, sin(2*M_PI*SystemFrequency*(p+1)*(j+k)*(1/RelaySamplingRate)+M_PI/2+2*M_PI/N));
                    col++;
                    gsl_matrix_set (A, j, col, cos(2*M_PI*SystemFrequency*(p+1)*(j+k)*(1/RelaySamplingRate)+M_PI/2+2*M_PI/N));
                    col++;}
                }

            gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, A, 0.0, AtA);
            int s;
            gsl_matrix_memcpy(tmp, AtA);
            gsl_linalg_LU_decomp(tmp, p, &s);
            gsl_linalg_LU_invert(tmp, p, inv_AtA);
            gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, inv_AtA, A, 0.0, PA);
            gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, PA, jan, 0.0, vX);

            recomp=2*Harmonic;
            imcomp=2*Harmonic+1;
            double retemp=M_SQRT1_2*gsl_matrix_get (vX, recomp, 0);
            double imtemp=M_SQRT1_2*gsl_matrix_get (vX, imcomp, 0);
            Sum= Complex(retemp,imtemp);
            Output->insert_Value(Sum);

			//cout<<abs(Output->get_Value())<<" "<<arg(Output->get_Value())<<" "<<real(Output->get_Value())<<" "<<imag(Output->get_Value())<<endl;

            k++;
            if(k>=2*N){k=0;}
		}

PLL::PLL(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,Channel<analogic> *CH_FREQ,float Freq){
			FunctionalDescription_set("Phase-Locked Loop");
			Input=CH_IN;
			Output=CH_OUT;
			SystemFrequency=Freq;
			Frequency=CH_FREQ;
			j=0;
		}
bool PLL::Prepare(float Samp){
			RelaySamplingRate=Samp;
			N=(unsigned)(RelaySamplingRate/SystemFrequency);
			u1=100; u2=10000; u3=0.02; C=0.001;
            A_old=1; W_old=2*M_PI*SystemFrequency; Fi_old=0;
            deltat=1/RelaySamplingRate;
			cout<<"PLL N"<<N<<endl;
			return true;
		}
void PLL::Run(){
			Complex Fasor,i(0.0,1.0);
			Fasor= Complex(0.0,0.0);

			sample=Input->get_Value();
			sig_pll=A_old*sin(Fi_old);
			Erro=sample-sig_pll;
			A_new=A_old+2*deltat*u1*Erro*sin(Fi_old);
			W_new=W_old+2*deltat*u2*Erro*A_old*(cos(Fi_old)/(A_old+C));
			Fi_new=Fi_old+deltat*W_old+2*deltat*u2*u3*Erro*A_old*(cos(Fi_old)/(A_old+C));

			Fasor_teta=Fi_new-W_new*j*deltat;
			Freq_med=W_new/(2*M_PI);
			Fasor_mod=A_new*M_SQRT1_2;

            Fi_old=Fi_new;
			W_old=W_new;
			A_old=A_new;

			float recomp=Fasor_mod*cos(Fasor_teta-M_PI/2-2*M_PI/N);
			float imcomp=Fasor_mod*sin(Fasor_teta-M_PI/2-2*M_PI/N);
			Fasor= Complex(recomp,imcomp);
            Output->insert_Value(Fasor);
            Frequency->insert_Value(Freq_med);

            //cout<<abs(Output->get_Value())<<" "<<arg(Output->get_Value())<<" "<<real(Output->get_Value())<<" "<<imag(Output->get_Value())<<endl;

			j++;
			if(j>=2*N){j=0;}
		}

RWT::RWT(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq){
			FunctionalDescription_set("Recursive Wavelet Transform");
			Input=CH_IN;
			Output=CH_OUT;
			SystemFrequency=Freq;
			k=0;
		}
bool RWT::Prepare(float Samp){
			RelaySamplingRate=Samp;
			N=(unsigned)(RelaySamplingRate/SystemFrequency);
			deltat=1/RelaySamplingRate;
			w0=2*M_PI;
			f=SystemFrequency;
			sig=2*M_PI/(sqrt(3));
			alfa_par=Complex(sig,-w0);
			alfa_par*=deltat*(-f);
			alfa=exp(alfa_par);

			W_1= Complex(0.0,0.0);
			W_2= Complex(0.0,0.0);
			W_3= Complex(0.0,0.0);
			W_4= Complex(0.0,0.0);
			W_5= Complex(0.0,0.0);
			W_6= Complex(0.0,0.0);

			temp=sig*f*deltat;

			lam1=alfa;
			lam1*=pow(temp,3)/3-pow(temp,4)/6+pow(temp,5)/15;
            lam2=pow(alfa,2);
            lam2*=(2*pow(temp,3)/3)-(5*pow(temp,4)/3)+(26*pow(temp,5)/15);
            lam3=pow(alfa,3);
            lam3*=(-6*pow(temp,3)/3)+(22*pow(temp,5)/5);
            lam4=pow(alfa,4);
            lam4*=(2*pow(temp,3)/3)+(5*pow(temp,4)/3)+(26*pow(temp,5)/15);
            lam5=pow(alfa,5);
            lam5*=(pow(temp,3)/3)+(pow(temp,4)/6)+(pow(temp,5)/15);

            beta1=alfa;
            beta1*=-6;
            beta2=pow(alfa,2);
            beta2*=15;
            beta3=pow(alfa,3);
            beta3*=-20;
            beta4=pow(alfa,4);
            beta4*=15;
            beta5=pow(alfa,5);
            beta5*=-6;
            beta6=pow(alfa,6);

            n1=Complex(sig,float(2*M_PI*f/f-w0));

			//cout<<"deltat:"<<deltat<<" w0:"<<w0<<" f:"<<f<<" sig:"<<sig<<" alfa:"<<alfa<<" lam1:"<<lam1<<" lam2:"<<lam2<<" lam3:"<<lam3<<" lam4:"<<lam4<<" lam5:"<<lam5<<" beta1:"<<beta1<<" beta2:"<<beta2<<" beta3:"<<beta3<<" beta4:"<<beta4<<" beta5:"<<beta5<<" beta6:"<<beta6<<" n1:"<<n1<<endl;
			cout<<"Recursive Wavelet N"<<N<<endl;
			return true;
		}
void RWT::Run(){
			Complex Fasor,i(0.0,1.0);
			Fasor= Complex(0.0,0.0);

			temp2=(k)*deltat*(-f);

			I_arg1=(-(((analogic)pow(temp2,3))/n1)*exp(n1*(analogic)temp2))+((((analogic)3*pow(temp2,2))/(pow(n1,2)))*exp(n1*(analogic)temp2))-((((analogic)6*temp2)/(pow(n1,3)))*exp(n1*(analogic)temp2))-((((analogic)-6)/(pow(n1,4)))*((analogic)1-exp(n1*(analogic)temp2)));
			I_arg2=(-(((analogic)pow(temp2,4))/n1)*exp(n1*(analogic)temp2))+((((analogic)4*pow(temp2,3))/(pow(n1,2)))*exp(n1*(analogic)temp2))-((((analogic)12*pow(temp2,2))/(pow(n1,3)))*exp(n1*(analogic)temp2))-((((analogic)24*temp2)/(pow(n1,4)))*exp(n1*(analogic)temp2))+((((analogic)24)/pow(n1,5))*((analogic)1-exp(n1*(analogic)temp2)));
			I_arg3=(-(((analogic)pow(temp2,5))/n1)*exp(n1*(analogic)temp2))+((((analogic)5*pow(temp2,4))/(pow(n1,2)))*exp(n1*(analogic)temp2))-((((analogic)20*pow(temp2,3))/(pow(n1,3)))*exp(n1*(analogic)temp2))+((((analogic)60*pow(temp2,2))/(pow(n1,4)))*exp(n1*(analogic)temp2))-((((analogic)120*temp2)/(pow(n1,5)))*exp(n1*(analogic)temp2))-(((analogic)120/(pow(n1,6)))*((analogic)1-exp(n1*(analogic)temp2)));

			I=((analogic)1/(sqrt(f)))*(((analogic)pow(-sig,3)/3)*(I_arg1)-((analogic)pow(sig,4)/6)*(I_arg2)-((analogic)pow(sig,5)/15)*(I_arg3));

            W=((analogic)deltat*sqrt(f))*(lam1*(analogic)Input->get_Value()+lam2*(analogic)Input->get_Value(1)+lam3*(analogic)Input->get_Value(2)+lam4*(analogic)Input->get_Value(3)+lam5*(analogic)Input->get_Value(4))-beta1*W_1-beta2*W_2-beta3*W_3-beta4*W_4-beta5*W_5-beta6*W_6;

            Am=abs(W/I)/(1.5*M_SQRT2);
            Fm=arg(W/I)-2*M_PI*f*deltat+M_PI_2;

            W_6=W_5;
            W_5=W_4;
            W_4=W_3;
            W_3=W_2;
            W_2=W_1;
            W_1=W;

            float recomp=Am*cos(Fm-M_PI/2);
			float imcomp=Am*sin(Fm-M_PI/2);
			Fasor= Complex(recomp,imcomp);
            Output->insert_Value(Fasor);

            //cout<<abs(Output->get_Value())<<" "<<arg(Output->get_Value())<<" "<<real(Output->get_Value())<<" "<<imag(Output->get_Value())<<endl;

			k++;
			if(k>=5051)
				k=0;
		}

RLS::RLS(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq){
			FunctionalDescription_set("Recursive Least Squares");
			Input=CH_IN;
			Output=CH_OUT;
			SystemFrequency=Freq;
        }
bool RLS::Prepare(float Samp){
			RelaySamplingRate=Samp;
			N=(unsigned)(RelaySamplingRate/SystemFrequency);
            ForgetingFactor=0.9;
            Ampi=0; Teta=0; deltat=1/RelaySamplingRate;
            k=0;

			cout<<"Recursive LS N"<<N<<endl;
			return true;
		}
void RLS::Run(){
			Complex Fasor,i(0.0,1.0);
            Fasor= Complex(0.0,0.0);

            gsl_matrix *xoi = gsl_matrix_calloc(2, 1);
            gsl_matrix *aoi = gsl_matrix_calloc(2, 1);
            gsl_matrix *Roi = gsl_matrix_calloc(2, 2);

            gsl_matrix *temp1 = gsl_matrix_calloc(2, 1);
            gsl_matrix *temp2 = gsl_matrix_calloc(2, 2);
            gsl_matrix *temp3 = gsl_matrix_calloc(2, 2);
            gsl_matrix *temp4 = gsl_matrix_calloc(1, 2);
            gsl_matrix *temp5 = gsl_matrix_calloc(1, 1);
            double temp6=0;
            gsl_matrix *temp7 = gsl_matrix_calloc(2, 1);

            fii=2*M_PI*SystemFrequency*k*deltat+Teta;
            e=-(Input->get_Value())+Ampi*sqrt(2)*sin(fii);
            gsl_matrix_set (xoi, 0, 0, sin(fii));
            gsl_matrix_set (xoi, 1, 0, Ampi*sqrt(2)*cos(fii));

            for (unsigned i=0;i<2;i++){
                    gsl_matrix_set(aoi,i,0, aoi_acum[i][0]);
                    for (unsigned j=0;j<2;j++){
                       gsl_matrix_set(Roi,i,j, Roi_acum[i][j]);
            }}

            if (k==0){
                gsl_matrix_set_identity (Roi);
                gsl_matrix_set_zero (aoi);
            }

            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Roi, xoi, 0.0, temp1);
            gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, temp1, xoi, 0.0, temp2);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp2, Roi, 0.0, temp3);
            gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, xoi, Roi, 0.0, temp4);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp4, xoi, 0.0, temp5);
            gsl_matrix_add_constant (temp5, ForgetingFactor);
            temp6=gsl_matrix_get(temp5,0,0);
            gsl_matrix_scale (temp3, 1/temp6);
            gsl_matrix_sub (Roi, temp3);
            gsl_matrix_scale (Roi, 1/ForgetingFactor);

            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Roi,xoi,0.0,temp7);
            gsl_matrix_scale (temp7, e);
            gsl_matrix_sub (aoi, temp7);

            for (unsigned i=0;i<2;i++){
                    aoi_acum[i][0]=gsl_matrix_get(aoi,i,0);
                    for (unsigned j=0;j<2;j++){
                        Roi_acum[i][j]=gsl_matrix_get(Roi,i,j);
                    }}

            Ampi=abs(gsl_matrix_get(aoi, 0, 0));
            Ampi/=sqrt(2);
            Teta=gsl_matrix_get(aoi, 1, 0);

            double recomp=Ampi*cos(Teta-M_PI/2);
            double imcomp=Ampi*sin(Teta-M_PI/2);
            Fasor=Complex (recomp,imcomp);

            Output->insert_Value(Fasor);

			//cout<<abs(Output->get_Value())<<" "<<arg(Output->get_Value())<<" "<<real(Output->get_Value())<<" "<<imag(Output->get_Value())<<endl;

            k++;
            if (k>=100*N){k=0;}
		}

DWT::DWT(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT, unsigned N, unsigned Family){
            FunctionalDescription_set("Discrete Wavelet Transform");
            Input=CH_IN;
            Output=CH_OUT;
            FamilyName=Family;
            NSamples=N;
            //Method=WAVELET_RMS;
            }
void DWT::Leafs(unsigned char ad,unsigned Nad,unsigned Level){
            Nivel.push_back(Level);
            AD.push_back(ad);
            Pos.push_back(Nad);
            }
void DWT::Run(){
            Complex Fasor,i(0.0,1.0);
            Fasor= Complex(0.0,0.0);

            double tempR1,tempR2;
            float teta1,teta2,teta,Mag;
            vector<double> AS;
            vector<double> AR3;
            double *S=new double[NSamples];
            neworder=true;

            for(unsigned k=0;k<NSamples;k++){
                S[k]=Input->get_Value(NSamples-1-k);
            }

            Transformada_Wavelet(S, NSamples, neworder);
            Transformada_Wavelet(S, NSamples/2, neworder);

            gsl_vector *R1_view = gsl_vector_calloc(NSamples/4);
            gsl_vector *R2_view = gsl_vector_calloc(NSamples/4);
            gsl_vector *S_view = gsl_vector_calloc(NSamples/4);
            for(unsigned j=0;j<NSamples/4;j++){
                AS.push_back(S[j]);
                gsl_vector_set(R1_view,j,AR1[j]);
                gsl_vector_set(R2_view,j,AR2[j]);
                gsl_vector_set(S_view,j,AS[j]);
                }

            double mAS=gsl_blas_dnrm2(S_view);

            gsl_blas_ddot(R1_view,S_view,&tempR1);
            gsl_blas_ddot(R2_view,S_view,&tempR2);
            teta1=acos(tempR1/(mAR1*mAS));
            teta2=acos(tempR2/(mAR2*mAS));

            if (teta2>=M_PI_2){teta=2*M_PI-teta1;}
            else {teta=teta1;}

            double *R3=new double[NSamples];
            for(unsigned j=0;j<NSamples;j++){
                R3[j]=sin(2*M_PI*j/NSamples+teta);
            }
            Transformada_Wavelet(R3, NSamples, neworder);
            Transformada_Wavelet(R3, NSamples/2, neworder);

            gsl_vector *R3_view = gsl_vector_calloc(NSamples/4);
            for(unsigned j=0;j<NSamples/4;j++){
                AR3.push_back(R3[j]);
                gsl_vector_set(R3_view,j,AR3[j]);
            }
            double mAR3=gsl_blas_dnrm2(R3_view);

            Mag=(mAS/mAR3)*M_SQRT1_2;

            float recomp=Mag*cos(teta);
            float imcomp=Mag*sin(teta);
            Fasor= Complex(recomp,imcomp);

            Output->insert_Value(Fasor);

            //cout<<abs(Output->get_Value())<<" "<<arg(Output->get_Value())<<" "<<real(Output->get_Value())<<" "<<imag(Output->get_Value())<<endl;

            delete(S);
            }
bool DWT::Prepare(float SampFreq){
                    unsigned Fim,Ini;
                    bool Inedito,Ord;
                    if(!getFilterBank()){
                        cout<<"Wavelet Family not Found"<<endl;
                        return false;
                    }
                    for(unsigned k=0;k<AD.size();k++){
                        cout<<"Leaf #"<<k<<endl;
                        if(Nivel[k]>(log(NSamples)/log(2))){
                            cout<<"Invalid Level Solicitation"<<endl;
                            return false;
                        }
                        if(Pos[k]>pow(2,Nivel[k]-1)){
                            cout<<"Invalid Detail/Aproximmation Solicitation"<<endl;
                            return false;
                        }
                        Ai.push_back(NSamples*((2*Pos[k]-2)/pow(2,Nivel[k])));
                        Af.push_back(NSamples*((2*Pos[k]-1)/pow(2,Nivel[k]))-1);
                        Di.push_back(NSamples*((2*Pos[k]-1)/pow(2,Nivel[k])));
                        Df.push_back(NSamples*((2*Pos[k])/pow(2,Nivel[k]))-1);
                        //cout<<"Ai("<<Pos[k]<<","<<Nivel[k]<<")="<<Ai[k]<<" Af("<<Pos[k]<<","<<Nivel[k]<<")="<<Af[k]<<endl;
                        //cout<<"Di("<<Pos[k]<<","<<Nivel[k]<<")="<<Di[k]<<" Df("<<Pos[k]<<","<<Nivel[k]<<")="<<Df[k]<<endl;
                        unsigned P=Pos[k];
                        for(unsigned i=Nivel[k];i>1;i--){
                            Inedito=true;
                            if((P%2)!=0){
                                P=(P+1)/2;
                                //cout<<"A"<<P<<","<<i-1<<" Direto ["<<NSamples*(2*P-2)/pow(2,i-1)<<","<<NSamples*(2*P-1)/pow(2,i-1)-1<<"]"<<endl;
                                Ini=NSamples*(2*P-2)/pow(2,i-1);
                                Fim=NSamples*(2*P-1)/pow(2,i-1)-1;
                                for(unsigned j=0;j<Ordem.size();j++){
                                    if((PosIni[j]==Ini)&&(PosFim[j]==Fim)){
                                    Inedito=false;
                                    }
                                }
                                if(Inedito){
                                    PosIni.push_back(Ini);
                                    PosFim.push_back(Fim);
                                    Ordem.push_back(true);
                                }
                            }
                            else{
                                P=P/2;
                                //cout<<"D"<<P<<","<<i-1<<" Inverso ["<<NSamples*(2*P-1)/pow(2,i-1)<<","<<NSamples*(2*P)/pow(2,i-1)-1<<"]"<<endl;
                                Ini=NSamples*(2*P-1)/pow(2,i-1);
                                Fim=NSamples*(2*P)/pow(2,i-1)-1;

                                for(unsigned j=0;j<Ordem.size();j++){
                                    if((PosIni[j]==Ini)&&(PosFim[j]==Fim)){
                                        Inedito=false;
                                    }
                                }
                                if(Inedito){
                                    PosIni.push_back(Ini);
                                    PosFim.push_back(Fim);
                                    Ordem.push_back(false);
                                }
                            }

                        }
                    }
                    //cout<<"I"<<" Direto ["<<0<<","<<NSamples-1<<"]"<<endl;
                    PosIni.push_back(0);
                    PosFim.push_back(NSamples-1);
                    Ordem.push_back(true);

                    for(unsigned k=0;k<Ordem.size();k++)
                        //cout<<"["<<PosIni[k]<<","<<PosFim[k]<<"] "<<Ordem[k]<<endl;
                    //cout<<endl;

                    for(unsigned k=0;k<Ordem.size();k++)
                    for(unsigned j=k;j<Ordem.size();j++)
                        if((PosFim[j]-PosIni[j])<(PosFim[k]-PosIni[k])){
                            Ini=PosIni[k]; Fim=PosFim[k];Ord=Ordem[k];
                            PosIni[k]=PosIni[j]; PosFim[k]=PosFim[j]; Ordem[k]=Ordem[j];
                            PosIni[j]=Ini; PosFim[j]=Fim; Ordem[j]=Ord;
                        }
                    for(unsigned k=0;k<Ordem.size();k++)
                        cout<<"["<<PosIni[k]<<","<<PosFim[k]<<"] "<<Ordem[k]<<endl;

                    //Alteração Prepare()//
                    //vector<double> R1;
                    double *R1=new double[NSamples];
                    double *R2=new double[NSamples];
                    //cout<<"-----------------------------"<<endl;
                    //cout<<"R1:  R2:    N:"<<NSamples<<endl;
                    for(unsigned j=0;j<NSamples;j++){
                        R1[j]=sin(2*M_PI*j/NSamples);
                        R2[j]=sin(2*M_PI*j/NSamples+M_PI_2);
                        //cout<<" "<<R1[j]<<"    "<<R2[j]<<"       "<<j<<endl;
                    }
                    //cout<<"-----------------------------"<<endl;
                    neworder=true;
                    Transformada_Wavelet(R1, NSamples, neworder);
                    Transformada_Wavelet(R1, NSamples/2, neworder);
                    Transformada_Wavelet(R2, NSamples, neworder);
                    Transformada_Wavelet(R2, NSamples/2, neworder);

                    for(unsigned j=0;j<NSamples/4;j++){
                        AR1.push_back(R1[j]);
                        AR2.push_back(R2[j]);
                    }
                    gsl_vector *R1_view = gsl_vector_calloc(AR1.size());
                    gsl_vector *R2_view = gsl_vector_calloc(AR2.size());

                    //cout<<"AR1:   AR2:"<<endl;
                    for(unsigned j=0;j<AR1.size();j++){
                        gsl_vector_set(R1_view,j,AR1[j]);
                        gsl_vector_set(R2_view,j,AR2[j]);
                        //cout<<" "<<AR1[j]<<"     "<<AR2[j]<<endl;
                    }
                    //cout<<"-----------------------------"<<endl;
                    mAR1=gsl_blas_dnrm2(R1_view);
                    mAR2=gsl_blas_dnrm2(R2_view);

                    //cout<<"mAR1:"<<mAR1<<" mAR2:"<<mAR2<<endl;
                    //cout<<"-----------------------------"<<endl;

                    cout<<"Discrete Wavelet N"<<NSamples<<endl;

                    return true;
                }
bool DWT::getFilterBank(){
                switch(FamilyName){
                    case haar:
                        FilterBank.push_back(0.7071067);FilterBank.push_back(0.7071067);
                        break;
                    case daub4:
                        FilterBank.push_back(4.82962913144e-01);FilterBank.push_back(8.36516303737e-01);FilterBank.push_back(2.24143868042e-01);FilterBank.push_back(-1.29409522551e-01);
                        break;
                    case daub6:
                        FilterBank.push_back(3.32670552950e-01);FilterBank.push_back(8.06891509311e-01);FilterBank.push_back(4.59877502118e-01);FilterBank.push_back(-1.35011020010e-01);FilterBank.push_back(-8.54412738820e-02);FilterBank.push_back(3.52262918857e-02);
                        break;
                    case daub8:
                        FilterBank.push_back(2.30377813308e-01);FilterBank.push_back(7.14846570552e-01);FilterBank.push_back(6.30880767929e-01);FilterBank.push_back(-2.79837694168e-02);FilterBank.push_back(-1.87034811719e-01);FilterBank.push_back(3.08413818355e-02);FilterBank.push_back(3.28830116668e-02);FilterBank.push_back(-1.05974017850e-02);
                        break;
                    case daub10:
                        FilterBank.push_back(1.60102397974e-01);FilterBank.push_back(6.03829269797e-01);FilterBank.push_back(7.24308528437e-01);FilterBank.push_back(1.38428145901e-01);FilterBank.push_back(-2.42294887066e-01);FilterBank.push_back(-3.22448695846e-02);FilterBank.push_back(7.75714938400e-02);FilterBank.push_back(-6.24149021279e-03);FilterBank.push_back(-1.25807519990e-02);FilterBank.push_back(3.33572528547e-03);
                        break;
                    case daub12:
                        FilterBank.push_back(1.11540743350e-01);FilterBank.push_back(4.94623890398e-01);FilterBank.push_back(7.51133908021e-01);FilterBank.push_back(3.15250351709e-01);FilterBank.push_back(-2.26264693965e-01);FilterBank.push_back(-1.29766867567e-01);FilterBank.push_back(9.75016055873e-02);FilterBank.push_back(2.75228655303e-02);FilterBank.push_back(-3.15820393174e-02);FilterBank.push_back(5.53842201161e-04);FilterBank.push_back(4.77725751094e-03);FilterBank.push_back(-1.07730108530e-03);
                        break;
                    case daub14:
                        FilterBank.push_back(7.78520540850e-02);FilterBank.push_back(3.96539319481e-01);FilterBank.push_back(7.29132090846e-01);FilterBank.push_back(4.69782287405e-01);FilterBank.push_back(-1.43906003928e-01);FilterBank.push_back(-2.24036184993e-01);FilterBank.push_back(7.13092192668e-02);FilterBank.push_back(8.06126091510e-02);FilterBank.push_back(-3.80299369350e-02);FilterBank.push_back(-1.65745416306e-02);FilterBank.push_back(1.25509985560e-02);FilterBank.push_back(4.29577972921e-04);FilterBank.push_back(-1.80164070404e-03);FilterBank.push_back(3.53713799974e-04);
                        break;
                    case daub16:
                        FilterBank.push_back(5.44158422431e-02);FilterBank.push_back(3.12871590914e-01);FilterBank.push_back(6.75630736297e-01);FilterBank.push_back(5.85354683654e-01);FilterBank.push_back(-1.58291052563e-02);FilterBank.push_back(-2.84015542961e-01);FilterBank.push_back(4.72484573913e-04);FilterBank.push_back(1.28747426620e-01);FilterBank.push_back(-1.73693010018e-02);FilterBank.push_back(-4.40882539307e-02);FilterBank.push_back(1.39810279173e-02);FilterBank.push_back(8.74609404740e-03);FilterBank.push_back(-4.87035299345e-03);FilterBank.push_back(-3.91740373376e-04);FilterBank.push_back(6.75449406450e-04);FilterBank.push_back(-1.17476784124e-04);
                        break;
                    case daub18:
                        FilterBank.push_back(3.80779473638e-02);FilterBank.push_back(2.43834674612e-01);FilterBank.push_back(6.04823123690e-01);FilterBank.push_back(6.57288078051e-01);FilterBank.push_back(1.33197385825e-01);FilterBank.push_back(-2.93273783279e-01);FilterBank.push_back(-9.68407832229e-02);FilterBank.push_back(1.48540749338e-01);FilterBank.push_back(3.07256814793e-02);FilterBank.push_back(-6.76328290613e-02);FilterBank.push_back(2.50947114831e-04);FilterBank.push_back(2.23616621236e-02);FilterBank.push_back(-4.72320475775e-03);FilterBank.push_back(-4.28150368246e-03);FilterBank.push_back(1.84764688305e-03);FilterBank.push_back(2.30385763523e-04);FilterBank.push_back(-2.51963188942e-04);FilterBank.push_back(3.93473203162e-05);
                        break;
                    case daub20:
                        FilterBank.push_back(2.66700579005e-02);FilterBank.push_back(1.88176800077e-01);FilterBank.push_back(5.27201188931e-01);FilterBank.push_back(6.88459039453e-01);FilterBank.push_back(2.81172343660e-01);FilterBank.push_back(-2.49846424327e-01);FilterBank.push_back(-1.95946274377e-01);FilterBank.push_back(1.27369340335e-01);FilterBank.push_back(9.30573646035e-02);FilterBank.push_back(-7.13941471663e-02);FilterBank.push_back(-2.94575368218e-02);FilterBank.push_back(3.32126740593e-02);FilterBank.push_back(3.60655356695e-03);FilterBank.push_back(-1.07331754833e-02);FilterBank.push_back(1.39535174705e-03);FilterBank.push_back(1.99240529518e-03);FilterBank.push_back(-6.85856694959e-04);FilterBank.push_back(-1.16466855129e-04);FilterBank.push_back(9.35886703200e-05);FilterBank.push_back(-1.32642028945e-05);
                        break;
                    case daub22:
                        FilterBank.push_back(1.86942977614e-02);FilterBank.push_back(1.44067021150e-01);FilterBank.push_back(4.49899764356e-01);FilterBank.push_back(6.85686774916e-01);FilterBank.push_back(4.11964368947e-01);FilterBank.push_back(-1.62275245027e-01);FilterBank.push_back(-2.74230846817e-01);FilterBank.push_back(6.60435881966e-02);FilterBank.push_back(1.49812012466e-01);FilterBank.push_back(-4.64799551166e-02);FilterBank.push_back(-6.64387856950e-02);FilterBank.push_back(3.13350902190e-02);FilterBank.push_back(2.08409043601e-02);FilterBank.push_back(-1.53648209062e-02);FilterBank.push_back(-3.34085887301e-03);FilterBank.push_back(4.92841765605e-03);FilterBank.push_back(-3.08592858815e-04);FilterBank.push_back(-8.93023250666e-04);FilterBank.push_back(2.49152523552e-04);FilterBank.push_back(5.44390746993e-05);FilterBank.push_back(-3.46349841869e-05);FilterBank.push_back(4.49427427723e-06);
                        break;
                    case daub24:
                        FilterBank.push_back(1.31122579572e-02);FilterBank.push_back(1.09566272821e-01);FilterBank.push_back(3.77355135214e-01);FilterBank.push_back(6.57198722579e-01);FilterBank.push_back(5.15886478427e-01);FilterBank.push_back(-4.47638856537e-02);FilterBank.push_back(-3.16178453752e-01);FilterBank.push_back(-2.37792572560e-02);FilterBank.push_back(1.82478605927e-01);FilterBank.push_back(5.35956967435e-03);FilterBank.push_back(-9.64321200965e-02);FilterBank.push_back(1.08491302558e-02);FilterBank.push_back(4.15462774950e-02);FilterBank.push_back(-1.22186490697e-02);FilterBank.push_back(-1.28408251983e-02);FilterBank.push_back(6.71149900879e-03);FilterBank.push_back(2.24860724099e-03);FilterBank.push_back(-2.17950361862e-03);FilterBank.push_back(6.54512821250e-06);FilterBank.push_back(3.88653062820e-04);FilterBank.push_back(-8.85041092082e-05);FilterBank.push_back(-2.42415457570e-05);FilterBank.push_back(1.27769522193e-05);FilterBank.push_back(-1.52907175806e-06);
                        break;
                    case daub26:
                        FilterBank.push_back(9.20213353896e-03);FilterBank.push_back(8.28612438729e-02);FilterBank.push_back(3.11996322160e-01);FilterBank.push_back(6.11055851158e-01);FilterBank.push_back(5.88889570431e-01);FilterBank.push_back(8.69857261796e-02);FilterBank.push_back(-3.14972907711e-01);FilterBank.push_back(-1.24576730750e-01);FilterBank.push_back(1.79476079429e-01);FilterBank.push_back(7.29489336567e-02);FilterBank.push_back(-1.05807618187e-01);FilterBank.push_back(-2.64884064753e-02);FilterBank.push_back(5.61394771002e-02);FilterBank.push_back(2.37997225405e-03);FilterBank.push_back(-2.38314207103e-02);FilterBank.push_back(3.92394144879e-03);FilterBank.push_back(7.25558940161e-03);FilterBank.push_back(-2.76191123465e-03);FilterBank.push_back(-1.31567391189e-03);FilterBank.push_back(9.32326130867e-04);FilterBank.push_back(4.92515251262e-05);FilterBank.push_back(-1.65128988556e-04);FilterBank.push_back(3.06785375793e-05);FilterBank.push_back(1.04419305714e-05);FilterBank.push_back(-4.70041647936e-06);FilterBank.push_back(5.22003509845e-07);
                        break;
                    case daub28:
                        FilterBank.push_back(6.46115346008e-03);FilterBank.push_back(6.23647588493e-02);FilterBank.push_back(2.54850267792e-01);FilterBank.push_back(5.54305617940e-01);FilterBank.push_back(6.31187849104e-01);FilterBank.push_back(2.18670687758e-01);FilterBank.push_back(-2.71688552278e-01);FilterBank.push_back(-2.18033529993e-01);FilterBank.push_back(1.38395213864e-01);FilterBank.push_back(1.39989016584e-01);FilterBank.push_back(-8.67484115681e-02);FilterBank.push_back(-7.15489555040e-02);FilterBank.push_back(5.52371262592e-02);FilterBank.push_back(2.69814083079e-02);FilterBank.push_back(-3.01853515403e-02);FilterBank.push_back(-5.61504953035e-03);FilterBank.push_back(1.27894932663e-02);FilterBank.push_back(-7.46218989268e-04);FilterBank.push_back(-3.84963886802e-03);FilterBank.push_back(1.06169108560e-03);FilterBank.push_back(7.08021154235e-04);FilterBank.push_back(-3.86831947312e-04);FilterBank.push_back(-4.17772457703e-05);FilterBank.push_back(6.87550425269e-05);FilterBank.push_back(-1.03372091845e-05);FilterBank.push_back(-4.38970490178e-06);FilterBank.push_back(1.72499467536e-06);FilterBank.push_back(-1.78713996831e-07);
                        break;
                    case daub30:
                        FilterBank.push_back(4.53853736157e-03);FilterBank.push_back(4.67433948927e-02);FilterBank.push_back(2.06023863986e-01);FilterBank.push_back(4.92631771708e-01);FilterBank.push_back(6.45813140357e-01);FilterBank.push_back(3.39002535454e-01);FilterBank.push_back(-1.93204139609e-01);FilterBank.push_back(-2.88882596566e-01);FilterBank.push_back(6.52829528487e-02);FilterBank.push_back(1.90146714007e-01);FilterBank.push_back(-3.96661765557e-02);FilterBank.push_back(-1.11120936037e-01);FilterBank.push_back(3.38771439235e-02);FilterBank.push_back(5.47805505845e-02);FilterBank.push_back(-2.57670073284e-02);FilterBank.push_back(-2.08100501696e-02);FilterBank.push_back(1.50839180278e-02);FilterBank.push_back(5.10100036040e-03);FilterBank.push_back(-6.48773456031e-03);FilterBank.push_back(-2.41756490761e-04);FilterBank.push_back(1.94332398038e-03);FilterBank.push_back(-3.73482354137e-04);FilterBank.push_back(-3.59565244362e-04);FilterBank.push_back(1.55896489920e-04);FilterBank.push_back(2.57926991553e-05);FilterBank.push_back(-2.81332962660e-05);FilterBank.push_back(3.36298718173e-06);FilterBank.push_back(1.81127040794e-06);FilterBank.push_back(-6.31688232588e-07);FilterBank.push_back(6.13335991330e-08);
                        break;
                    case daub32:
                        FilterBank.push_back(3.18922092534e-03);FilterBank.push_back(3.49077143236e-02);FilterBank.push_back(1.65064283488e-01);FilterBank.push_back(4.30312722846e-01);FilterBank.push_back(6.37356332083e-01);FilterBank.push_back(4.40290256886e-01);FilterBank.push_back(-8.97510894024e-02);FilterBank.push_back(-3.27063310527e-01);FilterBank.push_back(-2.79182081330e-02);FilterBank.push_back(2.11190693947e-01);FilterBank.push_back(2.73402637527e-02);FilterBank.push_back(-1.32388305563e-01);FilterBank.push_back(-6.23972275247e-03);FilterBank.push_back(7.59242360442e-02);FilterBank.push_back(-7.58897436885e-03);FilterBank.push_back(-3.68883976917e-02);FilterBank.push_back(1.02976596409e-02);FilterBank.push_back(1.39937688598e-02);FilterBank.push_back(-6.99001456341e-03);FilterBank.push_back(-3.64427962149e-03);FilterBank.push_back(3.12802338120e-03);FilterBank.push_back(4.07896980849e-04);FilterBank.push_back(-9.41021749359e-04);FilterBank.push_back(1.14241520038e-04);FilterBank.push_back(1.74787245225e-04);FilterBank.push_back(-6.10359662141e-05);FilterBank.push_back(-1.39456689882e-05);FilterBank.push_back(1.13366086612e-05);FilterBank.push_back(-1.04357134231e-06);FilterBank.push_back(-7.36365678545e-07);FilterBank.push_back(2.30878408685e-07);FilterBank.push_back(-2.10933963010e-08);
                        break;
                    case daub34:
                        FilterBank.push_back(2.24180700103e-03);FilterBank.push_back(2.59853937036e-02);FilterBank.push_back(1.31214903307e-01);FilterBank.push_back(3.70350724152e-01);FilterBank.push_back(6.10996615684e-01);FilterBank.push_back(5.18315764056e-01);FilterBank.push_back(2.73149704032e-02);FilterBank.push_back(-3.28320748363e-01);FilterBank.push_back(-1.26599752215e-01);FilterBank.push_back(1.97310589565e-01);FilterBank.push_back(1.01135489177e-01);FilterBank.push_back(-1.26815691778e-01);FilterBank.push_back(-5.70914196316e-02);FilterBank.push_back(8.11059866541e-02);FilterBank.push_back(2.23123361781e-02);FilterBank.push_back(-4.69224383892e-02);FilterBank.push_back(-3.27095553581e-03);FilterBank.push_back(2.27336765839e-02);FilterBank.push_back(-3.04298998135e-03);FilterBank.push_back(-8.60292152032e-03);FilterBank.push_back(2.96799669152e-03);FilterBank.push_back(2.30120524215e-03);FilterBank.push_back(-1.43684530480e-03);FilterBank.push_back(-3.28132519409e-04);FilterBank.push_back(4.39465427768e-04);FilterBank.push_back(-2.56101095665e-05);FilterBank.push_back(-8.20480320245e-05);FilterBank.push_back(2.31868137987e-05);FilterBank.push_back(6.99060098507e-06);FilterBank.push_back(-4.50594247722e-06);FilterBank.push_back(3.01654960999e-07);FilterBank.push_back(2.95770093331e-07);FilterBank.push_back(-8.42394844600e-08);FilterBank.push_back(7.26749296856e-09);
                        break;
                    case daub36:
                        FilterBank.push_back(1.57631021844e-03);FilterBank.push_back(1.92885317241e-02);FilterBank.push_back(1.03588465822e-01);FilterBank.push_back(3.14678941337e-01);FilterBank.push_back(5.71826807766e-01);FilterBank.push_back(5.71801654888e-01);FilterBank.push_back(1.47223111969e-01);FilterBank.push_back(-2.93654040736e-01);FilterBank.push_back(-2.16480934005e-01);FilterBank.push_back(1.49533975565e-01);FilterBank.push_back(1.67081312763e-01);FilterBank.push_back(-9.23318841508e-02);FilterBank.push_back(-1.06752246659e-01);FilterBank.push_back(6.48872162119e-02);FilterBank.push_back(5.70512477385e-02);FilterBank.push_back(-4.45261419029e-02);FilterBank.push_back(-2.37332103958e-02);FilterBank.push_back(2.66707059264e-02);FilterBank.push_back(6.26216795430e-03);FilterBank.push_back(-1.30514809466e-02);FilterBank.push_back(1.18630033858e-04);FilterBank.push_back(4.94334360546e-03);FilterBank.push_back(-1.11873266699e-03);FilterBank.push_back(-1.34059629833e-03);FilterBank.push_back(6.28465682965e-04);FilterBank.push_back(2.13581561910e-04);FilterBank.push_back(-1.98648552311e-04);FilterBank.push_back(-1.53591712353e-07);FilterBank.push_back(3.74123788074e-05);FilterBank.push_back(-8.52060253744e-06);FilterBank.push_back(-3.33263447888e-06);FilterBank.push_back(1.76871298362e-06);FilterBank.push_back(-7.69163268988e-08);FilterBank.push_back(-1.17609876702e-07);FilterBank.push_back(3.06883586304e-08);FilterBank.push_back(-2.50793445494e-09);
                        break;
                    case daub38:
                        FilterBank.push_back(1.10866976318e-03);FilterBank.push_back(1.42810984507e-02);FilterBank.push_back(8.12781132654e-02);FilterBank.push_back(2.64388431740e-01);FilterBank.push_back(5.24436377464e-01);FilterBank.push_back(6.01704549127e-01);FilterBank.push_back(2.60894952651e-01);FilterBank.push_back(-2.28091394215e-01);FilterBank.push_back(-2.85838631755e-01);FilterBank.push_back(7.46522697081e-02);FilterBank.push_back(2.12349743306e-01);FilterBank.push_back(-3.35185419023e-02);FilterBank.push_back(-1.42785695038e-01);FilterBank.push_back(2.75843506256e-02);FilterBank.push_back(8.69067555558e-02);FilterBank.push_back(-2.65012362501e-02);FilterBank.push_back(-4.56742262772e-02);FilterBank.push_back(2.16237674095e-02);FilterBank.push_back(1.93755498891e-02);FilterBank.push_back(-1.39883886785e-02);FilterBank.push_back(-5.86692228101e-03);FilterBank.push_back(7.04074736710e-03);FilterBank.push_back(7.68954359257e-04);FilterBank.push_back(-2.68755180070e-03);FilterBank.push_back(3.41808653458e-04);FilterBank.push_back(7.35802520505e-04);FilterBank.push_back(-2.60676135678e-04);FilterBank.push_back(-1.24600791734e-04);FilterBank.push_back(8.71127046721e-05);FilterBank.push_back(5.10595048707e-06);FilterBank.push_back(-1.66401762971e-05);FilterBank.push_back(3.01096431629e-06);FilterBank.push_back(1.53193147669e-06);FilterBank.push_back(-6.86275565776e-07);FilterBank.push_back(1.44708829879e-08);FilterBank.push_back(4.63693777578e-08);FilterBank.push_back(-1.11640206703e-08);FilterBank.push_back(8.66684883899e-10);
                        break;
                    case daub40:
                        FilterBank.push_back(7.79953613666e-04);FilterBank.push_back(1.05493946249e-02);FilterBank.push_back(6.34237804590e-02);FilterBank.push_back(2.19942113551e-01);FilterBank.push_back(4.72696185310e-01);FilterBank.push_back(6.10493238938e-01);FilterBank.push_back(3.61502298739e-01);FilterBank.push_back(-1.39212088011e-01);FilterBank.push_back(-3.26786800434e-01);FilterBank.push_back(-1.67270883090e-02);FilterBank.push_back(2.28291050819e-01);FilterBank.push_back(3.98502464577e-02);FilterBank.push_back(-1.55458750707e-01);FilterBank.push_back(-2.47168273386e-02);FilterBank.push_back(1.02291719174e-01);FilterBank.push_back(5.63224685730e-03);FilterBank.push_back(-6.17228996246e-02);FilterBank.push_back(5.87468181181e-03);FilterBank.push_back(3.22942995307e-02);FilterBank.push_back(-8.78932492390e-03);FilterBank.push_back(-1.38105261371e-02);FilterBank.push_back(6.72162730225e-03);FilterBank.push_back(4.42054238704e-03);FilterBank.push_back(-3.58149425960e-03);FilterBank.push_back(-8.31562172822e-04);FilterBank.push_back(1.39255961932e-03);FilterBank.push_back(-5.34975984399e-05);FilterBank.push_back(-3.85104748699e-04);FilterBank.push_back(1.01532889736e-04);FilterBank.push_back(6.77428082837e-05);FilterBank.push_back(-3.71058618339e-05);FilterBank.push_back(-4.37614386218e-06);FilterBank.push_back(7.24124828767e-06);FilterBank.push_back(-1.01199401001e-06);FilterBank.push_back(-6.84707959700e-07);FilterBank.push_back(2.63392422627e-07);FilterBank.push_back(2.01432202355e-10);FilterBank.push_back(-1.81484324829e-08);FilterBank.push_back(4.05612705555e-09);FilterBank.push_back(-2.99883648961e-10);
                        break;
                    case daub42:
                        FilterBank.push_back(5.48822509852e-04);FilterBank.push_back(7.77663905235e-03);FilterBank.push_back(4.92477715381e-02);FilterBank.push_back(1.81359625440e-01);FilterBank.push_back(4.19687944939e-01);FilterBank.push_back(6.01506094935e-01);FilterBank.push_back(4.44590451927e-01);FilterBank.push_back(-3.57229196172e-02);FilterBank.push_back(-3.35664089530e-01);FilterBank.push_back(-1.12397071568e-01);FilterBank.push_back(2.11564527680e-01);FilterBank.push_back(1.15233298439e-01);FilterBank.push_back(-1.39940424932e-01);FilterBank.push_back(-8.17759429808e-02);FilterBank.push_back(9.66003903237e-02);FilterBank.push_back(4.57234057492e-02);FilterBank.push_back(-6.49775048937e-02);FilterBank.push_back(-1.86538592021e-02);FilterBank.push_back(3.97268354278e-02);FilterBank.push_back(3.35775639033e-03);FilterBank.push_back(-2.08920536779e-02);FilterBank.push_back(2.40347092080e-03);FilterBank.push_back(8.98882438197e-03);FilterBank.push_back(-2.89133434858e-03);FilterBank.push_back(-2.95837403893e-03);FilterBank.push_back(1.71660704063e-03);FilterBank.push_back(6.39418500512e-04);FilterBank.push_back(-6.90671117082e-04);FilterBank.push_back(-3.19640627768e-05);FilterBank.push_back(1.93664650416e-04);FilterBank.push_back(-3.63552025008e-05);FilterBank.push_back(-3.49966598498e-05);FilterBank.push_back(1.53548250927e-05);FilterBank.push_back(2.79033053981e-06);FilterBank.push_back(-3.09001716454e-06);FilterBank.push_back(3.16609544236e-07);FilterBank.push_back(2.99213663046e-07);FilterBank.push_back(-1.00040087903e-07);FilterBank.push_back(-2.25401497467e-09);FilterBank.push_back(7.05803354123e-09);FilterBank.push_back(-1.47195419765e-09);FilterBank.push_back(1.03880557102e-10);
                        break;
                    case daub44:
                        FilterBank.push_back(3.86263231491e-04);FilterBank.push_back(5.72185463133e-03);FilterBank.push_back(3.80699372364e-02);FilterBank.push_back(1.48367540890e-01);FilterBank.push_back(3.67728683446e-01);FilterBank.push_back(5.78432731009e-01);FilterBank.push_back(5.07901090622e-01);FilterBank.push_back(7.37245011836e-02);FilterBank.push_back(-3.12726580428e-01);FilterBank.push_back(-2.00568406104e-01);FilterBank.push_back(1.64093188106e-01);FilterBank.push_back(1.79973187992e-01);FilterBank.push_back(-9.71107984091e-02);FilterBank.push_back(-1.31768137686e-01);FilterBank.push_back(6.80763143927e-02);FilterBank.push_back(8.45573763668e-02);FilterBank.push_back(-5.13642542974e-02);FilterBank.push_back(-4.65308118275e-02);FilterBank.push_back(3.69708466206e-02);FilterBank.push_back(2.05867076275e-02);FilterBank.push_back(-2.34800013444e-02);FilterBank.push_back(-6.21378284936e-03);FilterBank.push_back(1.25647252183e-02);FilterBank.push_back(3.00137398507e-04);FilterBank.push_back(-5.45569198615e-03);FilterBank.push_back(1.04426073918e-03);FilterBank.push_back(1.82701049565e-03);FilterBank.push_back(-7.70690988123e-04);FilterBank.push_back(-4.23787399839e-04);FilterBank.push_back(3.28609414213e-04);FilterBank.push_back(4.34589990453e-05);FilterBank.push_back(-9.40522363481e-05);FilterBank.push_back(1.13743496621e-05);FilterBank.push_back(1.73737569575e-05);FilterBank.push_back(-6.16672931646e-06);FilterBank.push_back(-1.56517913199e-06);FilterBank.push_back(1.29518205731e-06);FilterBank.push_back(-8.77987987336e-08);FilterBank.push_back(-1.28333622875e-07);FilterBank.push_back(3.76122874933e-08);FilterBank.push_back(1.68017140492e-09);FilterBank.push_back(-2.72962314663e-09);FilterBank.push_back(5.33593882166e-10);FilterBank.push_back(-3.60211348433e-11);
                        break;
                    case daub46:
                        FilterBank.push_back(2.71904194128e-04);FilterBank.push_back(4.20274889318e-03);FilterBank.push_back(2.93100036578e-02);FilterBank.push_back(1.20515531783e-01);FilterBank.push_back(3.18450813852e-01);FilterBank.push_back(5.44931147873e-01);FilterBank.push_back(5.51018517241e-01);FilterBank.push_back(1.81392625363e-01);FilterBank.push_back(-2.61392148030e-01);FilterBank.push_back(-2.71402098607e-01);FilterBank.push_back(9.21254070824e-02);FilterBank.push_back(2.23573658242e-01);FilterBank.push_back(-3.30374470942e-02);FilterBank.push_back(-1.64011321531e-01);FilterBank.push_back(2.02830745756e-02);FilterBank.push_back(1.12297043618e-01);FilterBank.push_back(-2.11262123562e-02);FilterBank.push_back(-7.02073915749e-02);FilterBank.push_back(2.17658568344e-02);FilterBank.push_back(3.84953325225e-02);FilterBank.push_back(-1.85235136501e-02);FilterBank.push_back(-1.75371010030e-02);FilterBank.push_back(1.27519439315e-02);FilterBank.push_back(6.03184065002e-03);FilterBank.push_back(-7.07531927370e-03);FilterBank.push_back(-1.13486547335e-03);FilterBank.push_back(3.12287644981e-03);FilterBank.push_back(-2.46501400516e-04);FilterBank.push_back(-1.06123122888e-03);FilterBank.push_back(3.19420492709e-04);FilterBank.push_back(2.56762452007e-04);FilterBank.push_back(-1.50021850349e-04);FilterBank.push_back(-3.37889483412e-05);FilterBank.push_back(4.42607120310e-05);FilterBank.push_back(-2.63520788924e-06);FilterBank.push_back(-8.34787556785e-06);FilterBank.push_back(2.39756954684e-06);FilterBank.push_back(8.14757483477e-07);FilterBank.push_back(-5.33900540520e-07);FilterBank.push_back(1.85309178563e-08);FilterBank.push_back(5.41754917953e-08);FilterBank.push_back(-1.39993549543e-08);FilterBank.push_back(-9.47288590181e-10);FilterBank.push_back(1.05044645369e-09);FilterBank.push_back(-1.93240511131e-10);FilterBank.push_back(1.25020330235e-11);
                        break;
                    case daub48:
                        FilterBank.push_back(1.91435800947e-04);FilterBank.push_back(3.08208171490e-03);FilterBank.push_back(2.24823399497e-02);FilterBank.push_back(9.72622358336e-02);FilterBank.push_back(2.72908916067e-01);FilterBank.push_back(5.04371040839e-01);FilterBank.push_back(5.74939221095e-01);FilterBank.push_back(2.80985553233e-01);FilterBank.push_back(-1.87271406885e-01);FilterBank.push_back(-3.17943078999e-01);FilterBank.push_back(4.77661368434e-03);FilterBank.push_back(2.39237388780e-01);FilterBank.push_back(4.25287296414e-02);FilterBank.push_back(-1.71175351370e-01);FilterBank.push_back(-3.87771735779e-02);FilterBank.push_back(1.21016303469e-01);FilterBank.push_back(2.09801137091e-02);FilterBank.push_back(-8.21616542080e-02);FilterBank.push_back(-4.57843624181e-03);FilterBank.push_back(5.13016200399e-02);FilterBank.push_back(-4.94470942812e-03);FilterBank.push_back(-2.82131070949e-02);FilterBank.push_back(7.66172188164e-03);FilterBank.push_back(1.30499708710e-02);FilterBank.push_back(-6.29143537001e-03);FilterBank.push_back(-4.74656878632e-03);FilterBank.push_back(3.73604617828e-03);FilterBank.push_back(1.15376493683e-03);FilterBank.push_back(-1.69645681897e-03);FilterBank.push_back(-4.41618485614e-05);FilterBank.push_back(5.86127059318e-04);FilterBank.push_back(-1.18123323796e-04);FilterBank.push_back(-1.46007981776e-04);FilterBank.push_back(6.55938863930e-05);FilterBank.push_back(2.18324146046e-05);FilterBank.push_back(-2.02288829261e-05);FilterBank.push_back(1.34115775080e-08);FilterBank.push_back(3.90110033859e-06);FilterBank.push_back(-8.98025314393e-07);FilterBank.push_back(-4.03250775687e-07);FilterBank.push_back(2.16633965327e-07);FilterBank.push_back(-5.05764541979e-10);FilterBank.push_back(-2.25574038817e-08);FilterBank.push_back(5.15777678967e-09);FilterBank.push_back(4.74837582425e-10);FilterBank.push_back(-4.02465864458e-10);FilterBank.push_back(6.99180115763e-11);FilterBank.push_back(-4.34278250380e-12);
                        break;
                    case daub50:
                        FilterBank.push_back(1.34802979347e-04);FilterBank.push_back(2.25695959185e-03);FilterBank.push_back(1.71867412540e-02);FilterBank.push_back(7.80358628721e-02);FilterBank.push_back(2.31693507886e-01);FilterBank.push_back(4.59683415146e-01);FilterBank.push_back(5.81636896746e-01);FilterBank.push_back(3.67885074802e-01);FilterBank.push_back(-9.71746409646e-02);FilterBank.push_back(-3.36473079641e-01);FilterBank.push_back(-8.75876145876e-02);FilterBank.push_back(2.24537819745e-01);FilterBank.push_back(1.18155286719e-01);FilterBank.push_back(-1.50560213750e-01);FilterBank.push_back(-9.85086152899e-02);FilterBank.push_back(1.06633805018e-01);FilterBank.push_back(6.67521644940e-02);FilterBank.push_back(-7.70841110565e-02);FilterBank.push_back(-3.71739628611e-02);FilterBank.push_back(5.36179093987e-02);FilterBank.push_back(1.55426059291e-02);FilterBank.push_back(-3.40423204606e-02);FilterBank.push_back(-3.07983679484e-03);FilterBank.push_back(1.89228044766e-02);FilterBank.push_back(-1.98942578220e-03);FilterBank.push_back(-8.86070261804e-03);FilterBank.push_back(2.72693625873e-03);FilterBank.push_back(3.32270777397e-03);FilterBank.push_back(-1.84248429020e-03);FilterBank.push_back(-8.99977423746e-04);FilterBank.push_back(8.77258193674e-04);FilterBank.push_back(1.15321244046e-04);FilterBank.push_back(-3.09880099098e-04);FilterBank.push_back(3.54371452327e-05);FilterBank.push_back(7.90464000396e-05);FilterBank.push_back(-2.73304811996e-05);FilterBank.push_back(-1.27719529319e-05);FilterBank.push_back(8.99066139306e-06);FilterBank.push_back(5.23282770815e-07);FilterBank.push_back(-1.77920133265e-06);FilterBank.push_back(3.21203751886e-07);FilterBank.push_back(1.92280679014e-07);FilterBank.push_back(-8.65694173227e-08);FilterBank.push_back(-2.61159855611e-09);FilterBank.push_back(9.27922448008e-09);FilterBank.push_back(-1.88041575506e-09);FilterBank.push_back(-2.22847491022e-10);FilterBank.push_back(1.53590157016e-10);FilterBank.push_back(-2.52762516346e-11);FilterBank.push_back(1.50969208282e-12);
                        break;
                    case daub52:
                        FilterBank.push_back(9.49379575071e-05);FilterBank.push_back(1.65052023353e-03);FilterBank.push_back(1.30975542925e-02);FilterBank.push_back(6.22747440251e-02);FilterBank.push_back(1.95039438716e-01);FilterBank.push_back(4.13292962278e-01);FilterBank.push_back(5.73669043034e-01);FilterBank.push_back(4.39158311789e-01);FilterBank.push_back(1.77407678098e-03);FilterBank.push_back(-3.26384593691e-01);FilterBank.push_back(-1.74839961289e-01);FilterBank.push_back(1.81291832311e-01);FilterBank.push_back(1.82755409589e-01);FilterBank.push_back(-1.04323900285e-01);FilterBank.push_back(-1.47977193275e-01);FilterBank.push_back(6.98231861132e-02);FilterBank.push_back(1.06482405249e-01);FilterBank.push_back(-5.34485616814e-02);FilterBank.push_back(-6.86547596040e-02);FilterBank.push_back(4.22321857963e-02);FilterBank.push_back(3.85357159711e-02);FilterBank.push_back(-3.13781103630e-02);FilterBank.push_back(-1.77609035683e-02);FilterBank.push_back(2.07349201799e-02);FilterBank.push_back(5.82958055531e-03);FilterBank.push_back(-1.17854979061e-02);FilterBank.push_back(-5.28738399262e-04);FilterBank.push_back(5.60194723942e-03);FilterBank.push_back(-9.39058250473e-04);FilterBank.push_back(-2.14553028156e-03);FilterBank.push_back(8.38348805654e-04);FilterBank.push_back(6.16138220457e-04);FilterBank.push_back(-4.31955707426e-04);FilterBank.push_back(-1.06057474828e-04);FilterBank.push_back(1.57479523860e-04);FilterBank.push_back(-5.27779549303e-06);FilterBank.push_back(-4.10967399639e-05);FilterBank.push_back(1.07422154087e-05);FilterBank.push_back(7.00007868296e-06);FilterBank.push_back(-3.88740016185e-06);FilterBank.push_back(-4.65046322064e-07);FilterBank.push_back(7.93921063370e-07);FilterBank.push_back(-1.07900423757e-07);FilterBank.push_back(-8.90446637016e-08);FilterBank.push_back(3.40779562129e-08);FilterBank.push_back(2.16932825985e-09);FilterBank.push_back(-3.77601047853e-09);FilterBank.push_back(6.78004724582e-10);FilterBank.push_back(1.00230319104e-10);FilterBank.push_back(-5.84040818534e-11);FilterBank.push_back(9.13051001637e-12);FilterBank.push_back(-5.25187122424e-13);
                        break;
                    case daub54:
                        FilterBank.push_back(6.68713138543e-05);FilterBank.push_back(1.20553123167e-03);FilterBank.push_back(9.95258878087e-03);FilterBank.push_back(4.94525999829e-02);FilterBank.push_back(1.62922027502e-01);FilterBank.push_back(3.67110214125e-01);FilterBank.push_back(5.53849860990e-01);FilterBank.push_back(4.93406122677e-01);FilterBank.push_back(1.02840855061e-01);FilterBank.push_back(-2.89716803314e-01);FilterBank.push_back(-2.48264581903e-01);FilterBank.push_back(1.14823019517e-01);FilterBank.push_back(2.27273288414e-01);FilterBank.push_back(-3.87864186318e-02);FilterBank.push_back(-1.78031740959e-01);FilterBank.push_back(1.57993974602e-02);FilterBank.push_back(1.31197971717e-01);FilterBank.push_back(-1.40627515558e-02);FilterBank.push_back(-9.10229065295e-02);FilterBank.push_back(1.73110182654e-02);FilterBank.push_back(5.79694057347e-02);FilterBank.push_back(-1.85124935619e-02);FilterBank.push_back(-3.27390666310e-02);FilterBank.push_back(1.61469669223e-02);FilterBank.push_back(1.56655956489e-02);FilterBank.push_back(-1.15771864589e-02);FilterBank.push_back(-5.86209634546e-03);FilterBank.push_back(6.85663560968e-03);FilterBank.push_back(1.34262687730e-03);FilterBank.push_back(-3.33285446952e-03);FilterBank.push_back(1.45752962593e-04);FilterBank.push_back(1.30117745024e-03);FilterBank.push_back(-3.41835122691e-04);FilterBank.push_back(-3.87901857410e-04);FilterBank.push_back(2.01971987969e-04);FilterBank.push_back(7.66005838706e-05);FilterBank.push_back(-7.71114551779e-05);FilterBank.push_back(-3.51748361490e-06);FilterBank.push_back(2.06344264773e-05);FilterBank.push_back(-3.90116407063e-06);FilterBank.push_back(-3.65750090818e-06);FilterBank.push_back(1.63436962472e-06);FilterBank.push_back(3.05088068625e-07);FilterBank.push_back(-3.47246814739e-07);FilterBank.push_back(3.28655896805e-08);FilterBank.push_back(4.02625505286e-08);FilterBank.push_back(-1.32133227399e-08);FilterBank.push_back(-1.30946560685e-09);FilterBank.push_back(1.52161498477e-09);FilterBank.push_back(-2.41552692801e-10);FilterBank.push_back(-4.37498622429e-11);FilterBank.push_back(2.21366208806e-11);FilterBank.push_back(-3.29579012247e-12);FilterBank.push_back(1.82818835288e-13);
                        break;
                    case daub56:
                        FilterBank.push_back(4.71080777501e-05);FilterBank.push_back(8.79498515984e-04);FilterBank.push_back(7.54265037764e-03);FilterBank.push_back(3.90926081154e-02);FilterBank.push_back(1.35137914253e-01);FilterBank.push_back(3.22563361285e-01);FilterBank.push_back(5.24998231630e-01);FilterBank.push_back(5.30516293441e-01);FilterBank.push_back(2.00176144045e-01);FilterBank.push_back(-2.30498954047e-01);FilterBank.push_back(-3.01327809532e-01);FilterBank.push_back(3.28578791633e-02);FilterBank.push_back(2.45808151373e-01);FilterBank.push_back(3.69068853157e-02);FilterBank.push_back(-1.82877330732e-01);FilterBank.push_back(-4.68382337445e-02);FilterBank.push_back(1.34627567910e-01);FilterBank.push_back(3.44786312750e-02);FilterBank.push_back(-9.76853558056e-02);FilterBank.push_back(-1.73419228313e-02);FilterBank.push_back(6.77478955019e-02);FilterBank.push_back(3.44801895554e-03);FilterBank.push_back(-4.33333686160e-02);FilterBank.push_back(4.43173291006e-03);FilterBank.push_back(2.46880600101e-02);FilterBank.push_back(-6.81554976455e-03);FilterBank.push_back(-1.20635919682e-02);FilterBank.push_back(5.83881662774e-03);FilterBank.push_back(4.78486311245e-03);FilterBank.push_back(-3.72546124707e-03);FilterBank.push_back(-1.36037384563e-03);FilterBank.push_back(1.87599866820e-03);FilterBank.push_back(1.41567239314e-04);FilterBank.push_back(-7.48674955911e-04);FilterBank.push_back(1.15465606365e-04);FilterBank.push_back(2.29579098223e-04);FilterBank.push_back(-8.90390149004e-05);FilterBank.push_back(-4.90771341619e-05);FilterBank.push_back(3.64140121105e-05);FilterBank.push_back(4.63866498139e-06);FilterBank.push_back(-1.00432604133e-05);FilterBank.push_back(1.24790031757e-06);FilterBank.push_back(1.84036373451e-06);FilterBank.push_back(-6.67021547995e-07);FilterBank.push_back(-1.75746117320e-07);FilterBank.push_back(1.49066001353e-07);FilterBank.push_back(-8.26238731562e-09);FilterBank.push_back(-1.78413869087e-08);FilterBank.push_back(5.04404705638e-09);FilterBank.push_back(6.94454032894e-10);FilterBank.push_back(-6.07704124722e-10);FilterBank.push_back(8.49222001105e-11);FilterBank.push_back(1.86736726378e-11);FilterBank.push_back(-8.36549047125e-12);FilterBank.push_back(1.18885053340e-12);FilterBank.push_back(-6.36777235471e-14);
                        break;
                    case daub58:
                        FilterBank.push_back(3.31896627984e-05);FilterBank.push_back(6.40951680304e-04);FilterBank.push_back(5.70212651777e-03);FilterBank.push_back(3.07735802214e-02);FilterBank.push_back(1.11370116951e-01);FilterBank.push_back(2.80653455970e-01);FilterBank.push_back(4.89758804762e-01);FilterBank.push_back(5.51374432758e-01);FilterBank.push_back(2.89105238335e-01);FilterBank.push_back(-1.54028734459e-01);FilterBank.push_back(-3.30040948917e-01);FilterBank.push_back(-5.57068000729e-02);FilterBank.push_back(2.36105236153e-01);FilterBank.push_back(1.12419174873e-01);FilterBank.push_back(-1.60877988594e-01);FilterBank.push_back(-1.07845949938e-01);FilterBank.push_back(1.14472295893e-01);FilterBank.push_back(8.32207471624e-02);FilterBank.push_back(-8.51254926156e-02);FilterBank.push_back(-5.50274895253e-02);FilterBank.push_back(6.34791645842e-02);FilterBank.push_back(3.05315432727e-02);FilterBank.push_back(-4.51879812777e-02);FilterBank.push_back(-1.29171425542e-02);FilterBank.push_back(2.94704318717e-02);FilterBank.push_back(2.64832730767e-03);FilterBank.push_back(-1.70412245736e-02);FilterBank.push_back(1.73788033272e-03);FilterBank.push_back(8.46972549356e-03);FilterBank.push_back(-2.55080712778e-03);FilterBank.push_back(-3.47379898968e-03);FilterBank.push_back(1.87712092572e-03);FilterBank.push_back(1.08705394222e-03);FilterBank.push_back(-1.00077832708e-03);FilterBank.push_back(-2.00071136307e-04);FilterBank.push_back(4.11128345474e-04);FilterBank.push_back(-2.29201804121e-05);FilterBank.push_back(-1.29304484008e-04);FilterBank.push_back(3.64502606856e-05);FilterBank.push_back(2.91334475016e-05);FilterBank.push_back(-1.65732839530e-05);FilterBank.push_back(-3.59364480402e-06);FilterBank.push_back(4.75060924645e-06);FilterBank.push_back(-3.02905459205e-07);FilterBank.push_back(-8.97570175063e-07);FilterBank.push_back(2.63389838699e-07);FilterBank.push_back(9.38719741109e-08);FilterBank.push_back(-6.28615692201e-08);FilterBank.push_back(1.07659190661e-09);FilterBank.push_back(7.76897885477e-09);FilterBank.push_back(-1.89399538617e-09);FilterBank.push_back(-3.42680086326e-10);FilterBank.push_back(2.40709945350e-10);FilterBank.push_back(-2.94058925076e-11);FilterBank.push_back(-7.83250973362e-12);FilterBank.push_back(3.15276241337e-12);FilterBank.push_back(-4.28565487006e-13);FilterBank.push_back(2.21919131158e-14);
                        break;
                    case daub60:
                        FilterBank.push_back(2.33861617273e-05);FilterBank.push_back(4.66637950428e-04);FilterBank.push_back(4.30079716504e-03);FilterBank.push_back(2.41308326715e-02);FilterBank.push_back(9.12383040670e-02);FilterBank.push_back(2.42020670940e-01);FilterBank.push_back(4.50487821853e-01);FilterBank.push_back(5.57572232912e-01);FilterBank.push_back(3.66242683371e-01);FilterBank.push_back(-6.61836707759e-02);FilterBank.push_back(-3.32966975020e-01);FilterBank.push_back(-1.41968513330e-01);FilterBank.push_back(1.99462121580e-01);FilterBank.push_back(1.77829873244e-01);FilterBank.push_back(-1.14558219432e-01);FilterBank.push_back(-1.57236817959e-01);FilterBank.push_back(7.27786589703e-02);FilterBank.push_back(1.22747746045e-01);FilterBank.push_back(-5.38064654582e-02);FilterBank.push_back(-8.76586900363e-02);FilterBank.push_back(4.38016646714e-02);FilterBank.push_back(5.67123657447e-02);FilterBank.push_back(-3.56733974967e-02);FilterBank.push_back(-3.22637589193e-02);FilterBank.push_back(2.70786195952e-02);FilterBank.push_back(1.52879607698e-02);FilterBank.push_back(-1.83997438681e-02);FilterBank.push_back(-5.29685966613e-03);FilterBank.push_back(1.09156316583e-02);FilterBank.push_back(6.19671756497e-04);FilterBank.push_back(-5.53073014819e-03);FilterBank.push_back(8.43384586662e-04);FilterBank.push_back(2.32452009406e-03);FilterBank.push_back(-8.60927696811e-04);FilterBank.push_back(-7.67878250438e-04);FilterBank.push_back(5.05094823903e-04);FilterBank.push_back(1.72482584235e-04);FilterBank.push_back(-2.16171830116e-04);FilterBank.push_back(-8.54830546758e-06);FilterBank.push_back(6.98200837080e-05);FilterBank.push_back(-1.33971686329e-05);FilterBank.push_back(-1.63615247872e-05);FilterBank.push_back(7.25214553589e-06);FilterBank.push_back(2.32754909849e-06);FilterBank.push_back(-2.18726767699e-06);FilterBank.push_back(1.09947433852e-08);FilterBank.push_back(4.26166232601e-07);FilterBank.push_back(-1.00041468235e-07);FilterBank.push_back(-4.76437996513e-08);FilterBank.push_back(2.60544275497e-08);FilterBank.push_back(5.55339786139e-10);FilterBank.push_back(-3.33110568046e-09);FilterBank.push_back(6.98486269183e-10);FilterBank.push_back(1.61362297827e-10);FilterBank.push_back(-9.46138799727e-11);FilterBank.push_back(1.00010513139e-11);FilterBank.push_back(3.23942863853e-12);FilterBank.push_back(-1.18523759210e-12);FilterBank.push_back(1.54399757084e-13);FilterBank.push_back(-7.73794263095e-15);
                        break;
                    case daub62:
                        FilterBank.push_back(1.64801338645e-05);FilterBank.push_back(3.39412203776e-04);FilterBank.push_back(3.23688406862e-03);FilterBank.push_back(1.88536916129e-02);FilterBank.push_back(7.43360930116e-02);FilterBank.push_back(2.07012874485e-01);FilterBank.push_back(4.09192200037e-01);FilterBank.push_back(5.51139840914e-01);FilterBank.push_back(4.29468808206e-01);FilterBank.push_back(2.71692124973e-02);FilterBank.push_back(-3.10955118319e-01);FilterBank.push_back(-2.17978485523e-01);FilterBank.push_back(1.40178288765e-01);FilterBank.push_back(2.24966711473e-01);FilterBank.push_back(-4.99263491604e-02);FilterBank.push_back(-1.86962360895e-01);FilterBank.push_back(1.54369884294e-02);FilterBank.push_back(1.45089500931e-01);FilterBank.push_back(-8.13983227346e-03);FilterBank.push_back(-1.07612773323e-01);FilterBank.push_back(1.09412974523e-02);FilterBank.push_back(7.53536117432e-02);FilterBank.push_back(-1.48800266181e-02);FilterBank.push_back(-4.86190754648e-02);FilterBank.push_back(1.61541715659e-02);FilterBank.push_back(2.80476193667e-02);FilterBank.push_back(-1.42762752777e-02);FilterBank.push_back(-1.39005529392e-02);FilterBank.push_back(1.05176394873e-02);FilterBank.push_back(5.51616357331e-03);FilterBank.push_back(-6.52085237587e-03);FilterBank.push_back(-1.42826422321e-03);FilterBank.push_back(3.39306677671e-03);FilterBank.push_back(-6.39790110601e-05);FilterBank.push_back(-1.45904174198e-03);FilterBank.push_back(3.43139829690e-04);FilterBank.push_back(4.99881617563e-04);FilterBank.push_back(-2.39658346940e-04);FilterBank.push_back(-1.24341161725e-04);FilterBank.push_back(1.08958435041e-04);FilterBank.push_back(1.50133572744e-05);FilterBank.push_back(-3.63125515786e-05);FilterBank.push_back(4.03452023518e-06);FilterBank.push_back(8.79530134269e-06);FilterBank.push_back(-3.03514236589e-06);FilterBank.push_back(-1.36906023094e-06);FilterBank.push_back(9.81001542204e-07);FilterBank.push_back(5.32725065697e-08);FilterBank.push_back(-1.97592512917e-07);FilterBank.push_back(3.61682651733e-08);FilterBank.push_back(2.32830971382e-08);FilterBank.push_back(-1.06152960215e-08);FilterBank.push_back(-6.47431168795e-10);FilterBank.push_back(1.40856815102e-09);FilterBank.push_back(-2.52404395415e-10);FilterBank.push_back(-7.34893003248e-11);FilterBank.push_back(3.69210880887e-11);FilterBank.push_back(-3.32700896712e-12);FilterBank.push_back(-1.32433491724e-12);FilterBank.push_back(4.44546709629e-13);FilterBank.push_back(-5.55944205057e-14);FilterBank.push_back(2.69938287976e-15);
                        break;
                    case daub64:
                        FilterBank.push_back(1.16146330213e-05);FilterBank.push_back(2.46656690638e-04);FilterBank.push_back(2.43126191957e-03);FilterBank.push_back(1.46810463814e-02);FilterBank.push_back(6.02574991203e-02);FilterBank.push_back(1.75750783639e-01);FilterBank.push_back(3.67509628597e-01);FilterBank.push_back(5.34317919340e-01);FilterBank.push_back(4.77809163733e-01);FilterBank.push_back(1.20630538265e-01);FilterBank.push_back(-2.66698181476e-01);FilterBank.push_back(-2.77421581558e-01);FilterBank.push_back(6.47133548055e-02);FilterBank.push_back(2.48310642356e-01);FilterBank.push_back(2.46624448396e-02);FilterBank.push_back(-1.92102344708e-01);FilterBank.push_back(-4.89951171846e-02);FilterBank.push_back(1.45232079475e-01);FilterBank.push_back(4.44049081999e-02);FilterBank.push_back(-1.09456113116e-01);FilterBank.push_back(-2.96278725084e-02);FilterBank.push_back(8.08741406384e-02);FilterBank.push_back(1.41061515161e-02);FilterBank.push_back(-5.69263140624e-02);FilterBank.push_back(-2.38026446493e-03);FilterBank.push_back(3.70514579235e-02);FilterBank.push_back(-4.14590766082e-03);FilterBank.push_back(-2.16628228363e-02);FilterBank.push_back(6.16752731068e-03);FilterBank.push_back(1.10174007154e-02);FilterBank.push_back(-5.41156825727e-03);FilterBank.push_back(-4.64921675118e-03);FilterBank.push_back(3.62722464068e-03);FilterBank.push_back(1.46895510046e-03);FilterBank.push_back(-1.96474055582e-03);FilterBank.push_back(-2.21167872957e-04);FilterBank.push_back(8.67305851845e-04);FilterBank.push_back(-1.02453731060e-04);FilterBank.push_back(-3.05965442382e-04);FilterBank.push_back(1.05391546173e-04);FilterBank.push_back(8.10367832913e-05);FilterBank.push_back(-5.25980928268e-05);FilterBank.push_back(-1.29404577940e-05);FilterBank.push_back(1.82426840198e-05);FilterBank.push_back(-6.36178153226e-07);FilterBank.push_back(-4.55830957626e-06);FilterBank.push_back(1.20288903632e-06);FilterBank.push_back(7.56004762559e-07);FilterBank.push_back(-4.28597069315e-07);FilterBank.push_back(-5.00336186874e-08);FilterBank.push_back(8.96596631195e-08);FilterBank.push_back(-1.21992435948e-08);FilterBank.push_back(-1.10438302172e-08);FilterBank.push_back(4.25042231198e-09);FilterBank.push_back(4.38438779994e-10);FilterBank.push_back(-5.88109146263e-10);FilterBank.push_back(8.90472379622e-11);FilterBank.push_back(3.26327074133e-11);FilterBank.push_back(-1.43091876516e-11);FilterBank.push_back(1.07561065350e-12);FilterBank.push_back(5.36148222961e-13);FilterBank.push_back(-1.66380048943e-13);FilterBank.push_back(2.00071530381e-14);FilterBank.push_back(-9.42101913953e-16);
                        break;
                    case daub66:
                        FilterBank.push_back(8.18635831417e-06);FilterBank.push_back(1.79101615370e-04);FilterBank.push_back(1.82270943516e-03);FilterBank.push_back(1.13959433745e-02);FilterBank.push_back(4.86146665317e-02);FilterBank.push_back(1.48186313180e-01);FilterBank.push_back(3.26718130117e-01);FilterBank.push_back(5.09376172514e-01);FilterBank.push_back(5.11254770583e-01);FilterBank.push_back(2.09582350713e-01);FilterBank.push_back(-2.04202622398e-01);FilterBank.push_back(-3.15997410766e-01);FilterBank.push_back(-1.92783394369e-02);FilterBank.push_back(2.45420612119e-01);FilterBank.push_back(9.98515586803e-02);FilterBank.push_back(-1.71428099051e-01);FilterBank.push_back(-1.10844133116e-01);FilterBank.push_back(1.21967856403e-01);FilterBank.push_back(9.47880880506e-02);FilterBank.push_back(-9.11469683513e-02);FilterBank.push_back(-7.03024850540e-02);FilterBank.push_back(7.01911439409e-02);FilterBank.push_back(4.57345618938e-02);FilterBank.push_back(-5.34712513358e-02);FilterBank.push_back(-2.52485829774e-02);FilterBank.push_back(3.86870607602e-02);FilterBank.push_back(1.07032658200e-02);FilterBank.push_back(-2.57287617547e-02);FilterBank.push_back(-2.16775861735e-03);FilterBank.push_back(1.53169541158e-02);FilterBank.push_back(-1.59428878241e-03);FilterBank.push_back(-7.95354038705e-03);FilterBank.push_back(2.38906240816e-03);FilterBank.push_back(3.48080095340e-03);FilterBank.push_back(-1.86071821445e-03);FilterBank.push_back(-1.20430925760e-03);FilterBank.push_back(1.07438069635e-03);FilterBank.push_back(2.72730584733e-04);FilterBank.push_back(-4.90832900759e-04);FilterBank.push_back(4.39316625176e-06);FilterBank.push_back(1.78043189825e-04);FilterBank.push_back(-4.16043851627e-05);FilterBank.push_back(-4.92956442341e-05);FilterBank.push_back(2.42333539881e-05);FilterBank.push_back(9.07080575782e-06);FilterBank.push_back(-8.86612136675e-06);FilterBank.push_back(-3.60751610287e-07);FilterBank.push_back(2.28837127614e-06);FilterBank.push_back(-4.42692340795e-07);FilterBank.push_back(-3.98579129198e-07);FilterBank.push_back(1.82244333257e-07);FilterBank.push_back(3.37797270373e-08);FilterBank.push_back(-3.98783819851e-08);FilterBank.push_back(3.67286357683e-09);FilterBank.push_back(5.11121185734e-09);FilterBank.push_back(-1.67139267725e-09);FilterBank.push_back(-2.49640210524e-10);FilterBank.push_back(2.42683310230e-10);FilterBank.push_back(-3.04957445394e-11);FilterBank.push_back(-1.42023685988e-11);FilterBank.push_back(5.50941472076e-12);FilterBank.push_back(-3.34348121895e-13);FilterBank.push_back(-2.15248838683e-13);FilterBank.push_back(6.21474024717e-14);FilterBank.push_back(-7.19651054536e-15);FilterBank.push_back(3.28937367841e-16);
                        break;
                    case daub68:
                        FilterBank.push_back(5.77051063273e-06);FilterBank.push_back(1.29947620067e-04);FilterBank.push_back(1.36406139005e-03);FilterBank.push_back(8.81988940388e-03);FilterBank.push_back(3.90488413517e-02);FilterBank.push_back(1.24152482111e-01);FilterBank.push_back(2.87765059233e-01);FilterBank.push_back(4.78478746279e-01);FilterBank.push_back(5.30555099656e-01);FilterBank.push_back(2.90366329507e-01);FilterBank.push_back(-1.28246842174e-01);FilterBank.push_back(-3.31525301508e-01);FilterBank.push_back(-1.03891915515e-01);FilterBank.push_back(2.16907220187e-01);FilterBank.push_back(1.66601750412e-01);FilterBank.push_back(-1.27337358223e-01);FilterBank.push_back(-1.60924927177e-01);FilterBank.push_back(7.79918469379e-02);FilterBank.push_back(1.34125960271e-01);FilterBank.push_back(-5.44829680641e-02);FilterBank.push_back(-1.02947596992e-01);FilterBank.push_back(4.35760946496e-02);FilterBank.push_back(7.31852354367e-02);FilterBank.push_back(-3.70128384178e-02);FilterBank.push_back(-4.74385596452e-02);FilterBank.push_back(3.07397465739e-02);FilterBank.push_back(2.72283507563e-02);FilterBank.push_back(-2.36717379228e-02);FilterBank.push_back(-1.31439800166e-02);FilterBank.push_back(1.64093741998e-02);FilterBank.push_back(4.71364926099e-03);FilterBank.push_back(-1.00455067083e-02);FilterBank.push_back(-6.19474884515e-04);FilterBank.push_back(5.33495076875e-03);FilterBank.push_back(-7.69212797506e-04);FilterBank.push_back(-2.39945394353e-03);FilterBank.push_back(8.58995987436e-04);FilterBank.push_back(8.75199906407e-04);FilterBank.push_back(-5.52735576214e-04);FilterBank.push_back(-2.32673214023e-04);FilterBank.push_back(2.65077239755e-04);FilterBank.push_back(2.66005001845e-05);FilterBank.push_back(-9.91469777078e-05);FilterBank.push_back(1.35311722724e-05);FilterBank.push_back(2.84495141969e-05);FilterBank.push_back(-1.05765749425e-05);FilterBank.push_back(-5.71082651099e-06);FilterBank.push_back(4.16987175854e-06);FilterBank.push_back(4.97971810142e-07);FilterBank.push_back(-1.11630653481e-06);FilterBank.push_back(1.44819570833e-07);FilterBank.push_back(2.02599066666e-07);FilterBank.push_back(-7.52670174041e-08);FilterBank.push_back(-1.99034650153e-08);FilterBank.push_back(1.74042333293e-08);FilterBank.push_back(-8.66574426136e-10);FilterBank.push_back(-2.31650194699e-09);FilterBank.push_back(6.44637821032e-10);FilterBank.push_back(1.30041031860e-10);FilterBank.push_back(-9.90477453763e-11);FilterBank.push_back(1.00420873546e-11);FilterBank.push_back(6.08012535400e-12);FilterBank.push_back(-2.10787910891e-12);FilterBank.push_back(9.79945115821e-14);FilterBank.push_back(8.57919405179e-14);FilterBank.push_back(-2.31708370390e-14);FilterBank.push_back(2.58733838193e-15);FilterBank.push_back(-1.14894475448e-16);
                        break;
                    case daub70:
                        FilterBank.push_back(4.06793406114e-06);FilterBank.push_back(9.42146947557e-05);FilterBank.push_back(1.01912268037e-03);FilterBank.push_back(6.80729288431e-03);FilterBank.push_back(3.12362885114e-02);FilterBank.push_back(1.03404455861e-01);FilterBank.push_back(2.51307378994e-01);FilterBank.push_back(4.43592739224e-01);FilterBank.push_back(5.37008427509e-01);FilterBank.push_back(3.60345640518e-01);FilterBank.push_back(-4.38838818739e-02);FilterBank.push_back(-3.23822864912e-01);FilterBank.push_back(-1.81786976766e-01);FilterBank.push_back(1.66041357490e-01);FilterBank.push_back(2.17299289321e-01);FilterBank.push_back(-6.52628713106e-02);FilterBank.push_back(-1.91919589298e-01);FilterBank.push_back(1.93095446660e-02);FilterBank.push_back(1.55292480396e-01);FilterBank.push_back(-4.75268083411e-03);FilterBank.push_back(-1.20585522643e-01);FilterBank.push_back(4.73422917264e-03);FilterBank.push_back(8.99135475707e-02);FilterBank.push_back(-9.31855894990e-03);FilterBank.push_back(-6.33560374404e-02);FilterBank.push_back(1.32285495850e-02);FilterBank.push_back(4.12546930647e-02);FilterBank.push_back(-1.43668397842e-02);FilterBank.push_back(-2.41694978016e-02);FilterBank.push_back(1.27664567156e-02);FilterBank.push_back(1.22894360081e-02);FilterBank.push_back(-9.57779789923e-03);FilterBank.push_back(-5.08599164923e-03);FilterBank.push_back(6.13775458674e-03);FilterBank.push_back(1.42808879407e-03);FilterBank.push_back(-3.35764438092e-03);FilterBank.push_back(7.61596943517e-06);FilterBank.push_back(1.54963746970e-03);FilterBank.push_back(-3.34669216425e-04);FilterBank.push_back(-5.86481031899e-04);FilterBank.push_back(2.64832881996e-04);FilterBank.push_back(1.70001228366e-04);FilterBank.push_back(-1.36588307226e-04);FilterBank.push_back(-2.97699596284e-05);FilterBank.push_back(5.30414312291e-05);FilterBank.push_back(-2.43700152682e-06);FilterBank.push_back(-1.57244207727e-05);FilterBank.push_back(4.30804786171e-06);FilterBank.push_back(3.35334586287e-06);FilterBank.push_back(-1.89592961769e-06);FilterBank.push_back(-3.90393173328e-07);FilterBank.push_back(5.30236861690e-07);FilterBank.push_back(-3.70030837820e-08);FilterBank.push_back(-9.99039694453e-08);FilterBank.push_back(3.00818865071e-08);FilterBank.push_back(1.08490273378e-08);FilterBank.push_back(-7.45811655289e-09);FilterBank.push_back(5.89795131038e-11);FilterBank.push_back(1.03082334548e-09);FilterBank.push_back(-2.43354557375e-10);FilterBank.push_back(-6.40793825650e-11);FilterBank.push_back(4.00053662725e-11);FilterBank.push_back(-3.12563935710e-12);FilterBank.push_back(-2.56706547615e-12);FilterBank.push_back(8.01508853368e-13);FilterBank.push_back(-2.59795432889e-14);FilterBank.push_back(-3.39772085679e-14);FilterBank.push_back(8.62403743472e-15);FilterBank.push_back(-9.29801252932e-16);FilterBank.push_back(4.01462871233e-17);
                        break;
                    case daub72:
                        FilterBank.push_back(2.86792518275e-06);FilterBank.push_back(6.82602867854e-05);FilterBank.push_back(7.60215109966e-04);FilterBank.push_back(5.24029737740e-03);FilterBank.push_back(2.48905656448e-02);FilterBank.push_back(8.56520925952e-02);FilterBank.push_back(2.17756953097e-01);FilterBank.push_back(4.06433697708e-01);FilterBank.push_back(5.32266895260e-01);FilterBank.push_back(4.17875335600e-01);FilterBank.push_back(4.39751975293e-02);FilterBank.push_back(-2.94421039589e-01);FilterBank.push_back(-2.46807036978e-01);FilterBank.push_back(9.81142041631e-02);FilterBank.push_back(2.46537277608e-01);FilterBank.push_back(7.27851509579e-03);FilterBank.push_back(-1.99337205608e-01);FilterBank.push_back(-4.58614007463e-02);FilterBank.push_back(1.54106236627e-01);FilterBank.push_back(5.02761800735e-02);FilterBank.push_back(-1.18803754310e-01);FilterBank.push_back(-3.98808535755e-02);FilterBank.push_back(9.11567822580e-02);FilterBank.push_back(2.50387214495e-02);FilterBank.push_back(-6.82090166368e-02);FilterBank.push_back(-1.13191003168e-02);FilterBank.push_back(4.85130835478e-02);FilterBank.push_back(1.42497266176e-03);FilterBank.push_back(-3.19807206776e-02);FilterBank.push_back(3.98404019871e-03);FilterBank.push_back(1.90635947806e-02);FilterBank.push_back(-5.65781324505e-03);FilterBank.push_back(-9.99026347328e-03);FilterBank.push_back(5.02298910666e-03);FilterBank.push_back(4.41348483535e-03);FilterBank.push_back(-3.48454144540e-03);FilterBank.push_back(-1.50307406629e-03);FilterBank.push_back(1.99079377185e-03);FilterBank.push_back(2.77681279571e-04);FilterBank.push_back(-9.46340382326e-04);FilterBank.push_back(8.61456575899e-05);FilterBank.push_back(3.69350728496e-04);FilterBank.push_back(-1.15511889584e-04);FilterBank.push_back(-1.13189946808e-04);FilterBank.push_back(6.69474119693e-05);FilterBank.push_back(2.37510668366e-05);FilterBank.push_back(-2.73139082465e-05);FilterBank.push_back(-1.18347105998e-06);FilterBank.push_back(8.37221819816e-06);FilterBank.push_back(-1.58614578243e-06);FilterBank.push_back(-1.87081160285e-06);FilterBank.push_back(8.31142127970e-07);FilterBank.push_back(2.54842352255e-07);FilterBank.push_back(-2.45537765843e-07);FilterBank.push_back(2.75324907333e-09);FilterBank.push_back(4.79904346545e-08);FilterBank.push_back(-1.15609368881e-08);FilterBank.push_back(-5.61278434332e-09);FilterBank.push_back(3.13884169578e-09);FilterBank.push_back(1.09081555371e-10);FilterBank.push_back(-4.51254577856e-10);FilterBank.push_back(8.96241820385e-11);FilterBank.push_back(3.03742909811e-11);FilterBank.push_back(-1.59971668926e-11);FilterBank.push_back(8.87684628721e-13);FilterBank.push_back(1.07096935711e-12);FilterBank.push_back(-3.02928502697e-13);FilterBank.push_back(5.54226318263e-15);FilterBank.push_back(1.33807138629e-14);FilterBank.push_back(-3.20462854340e-15);FilterBank.push_back(3.33997198481e-16);FilterBank.push_back(-1.40327417537e-17);
                        break;
                    case daub74:
                        FilterBank.push_back(2.02206086249e-06);FilterBank.push_back(4.94234375062e-05);FilterBank.push_back(5.66241837706e-04);FilterBank.push_back(4.02414036825e-03);FilterBank.push_back(1.97622861538e-02);FilterBank.push_back(7.05848259771e-02);FilterBank.push_back(1.87326331862e-01);FilterBank.push_back(3.68440972400e-01);FilterBank.push_back(5.18167040855e-01);FilterBank.push_back(4.62207553661e-01);FilterBank.push_back(1.30878963233e-01);FilterBank.push_back(-2.46180429761e-01);FilterBank.push_back(-2.94375915262e-01);FilterBank.push_back(1.96715004523e-02);FilterBank.push_back(2.51523254360e-01);FilterBank.push_back(8.18060283872e-02);FilterBank.push_back(-1.81962291778e-01);FilterBank.push_back(-1.08451713823e-01);FilterBank.push_back(1.29929646959e-01);FilterBank.push_back(1.01780296838e-01);FilterBank.push_back(-9.66075406166e-02);FilterBank.push_back(-8.23302119065e-02);FilterBank.push_back(7.50476199483e-02);FilterBank.push_back(5.95674108715e-02);FilterBank.push_back(-5.92568156326e-02);FilterBank.push_back(-3.82538294793e-02);FilterBank.push_back(4.58079441512e-02);FilterBank.push_back(2.09728005925e-02);FilterBank.push_back(-3.35235840641e-02);FilterBank.push_back(-8.83349389041e-03);FilterBank.push_back(2.26186515445e-02);FilterBank.push_back(1.69047238348e-03);FilterBank.push_back(-1.37639819628e-02);FilterBank.push_back(1.51930577883e-03);FilterBank.push_back(7.38775745285e-03);FilterBank.push_back(-2.24805318700e-03);FilterBank.push_back(-3.39452327640e-03);FilterBank.push_back(1.81687134380e-03);FilterBank.push_back(1.26393425811e-03);FilterBank.push_back(-1.11148486531e-03);FilterBank.push_back(-3.28078847088e-04);FilterBank.push_back(5.49053277337e-04);FilterBank.push_back(1.53443902319e-05);FilterBank.push_back(-2.20894403245e-04);FilterBank.push_back(4.33672612594e-05);FilterBank.push_back(7.05513878206e-05);FilterBank.push_back(-3.09866292761e-05);FilterBank.push_back(-1.63916249616e-05);FilterBank.push_back(1.35432771841e-05);FilterBank.push_back(1.84994500311e-06);FilterBank.push_back(-4.30994155659e-06);FilterBank.push_back(4.85473139699e-07);FilterBank.push_back(1.00212139929e-06);FilterBank.push_back(-3.49494860344e-07);FilterBank.push_back(-1.50988538867e-07);FilterBank.push_back(1.10903123221e-07);FilterBank.push_back(5.35065751546e-09);FilterBank.push_back(-2.25219383672e-08);FilterBank.push_back(4.22448570636e-09);FilterBank.push_back(2.79397446595e-09);FilterBank.push_back(-1.29720500146e-09);FilterBank.push_back(-1.03141112909e-10);FilterBank.push_back(1.94616489408e-10);FilterBank.push_back(-3.20339824412e-11);FilterBank.push_back(-1.39841571553e-11);FilterBank.push_back(6.33495544097e-12);FilterBank.push_back(-2.09636319423e-13);FilterBank.push_back(-4.42161240987e-13);FilterBank.push_back(1.13805283092e-13);FilterBank.push_back(-4.51888960746e-16);FilterBank.push_back(-5.24302569188e-15);FilterBank.push_back(1.18901238750e-15);FilterBank.push_back(-1.19928033585e-16);FilterBank.push_back(4.90661506493e-18);
                        break;
                    case daub76:
                        FilterBank.push_back(1.42577664167e-06);FilterBank.push_back(3.57625199426e-05);FilterBank.push_back(4.21170266472e-04);FilterBank.push_back(3.08308811925e-03);FilterBank.push_back(1.56372493475e-02);FilterBank.push_back(5.78899436128e-02);FilterBank.push_back(1.60071993564e-01);FilterBank.push_back(3.30775781411e-01);FilterBank.push_back(4.96591175311e-01);FilterBank.push_back(4.93356078517e-01);FilterBank.push_back(2.13050571355e-01);FilterBank.push_back(-1.82867667708e-01);FilterBank.push_back(-3.21675637808e-01);FilterBank.push_back(-6.22665060478e-02);FilterBank.push_back(2.32125963835e-01);FilterBank.push_back(1.49985119618e-01);FilterBank.push_back(-1.41795685973e-01);FilterBank.push_back(-1.59912565158e-01);FilterBank.push_back(8.56381215561e-02);FilterBank.push_back(1.41414734073e-01);FilterBank.push_back(-5.65864586307e-02);FilterBank.push_back(-1.14731170710e-01);FilterBank.push_back(4.30958954330e-02);FilterBank.push_back(8.72043982620e-02);FilterBank.push_back(-3.66051034028e-02);FilterBank.push_back(-6.17662087084e-02);FilterBank.push_back(3.19898775315e-02);FilterBank.push_back(4.00549811051e-02);FilterBank.push_back(-2.68914938808e-02);FilterBank.push_back(-2.31141340205e-02);FilterBank.push_back(2.09046452556e-02);FilterBank.push_back(1.12904972786e-02);FilterBank.push_back(-1.47018820653e-02);FilterBank.push_back(-4.13130665603e-03);FilterBank.push_back(9.21478503219e-03);FilterBank.push_back(5.62571574840e-04);FilterBank.push_back(-5.07131450921e-03);FilterBank.push_back(7.16982182106e-04);FilterBank.push_back(2.40069778189e-03);FilterBank.push_back(-8.44862666553e-04);FilterBank.push_back(-9.42461407722e-04);FilterBank.push_back(5.81075975053e-04);FilterBank.push_back(2.81763925038e-04);FilterBank.push_back(-3.03102046072e-04);FilterBank.push_back(-4.55568269666e-05);FilterBank.push_back(1.26204335016e-04);FilterBank.push_back(-1.15540910383e-05);FilterBank.push_back(-4.17514164854e-05);FilterBank.push_back(1.33417614992e-05);FilterBank.push_back(1.03735918404e-05);FilterBank.push_back(-6.45673042846e-06);FilterBank.push_back(-1.55084435011e-06);FilterBank.push_back(2.14996026993e-06);FilterBank.push_back(-8.48708758607e-08);FilterBank.push_back(-5.18773373887e-07);FilterBank.push_back(1.39637754550e-07);FilterBank.push_back(8.40035104689e-08);FilterBank.push_back(-4.88475793745e-08);FilterBank.push_back(-5.42427480028e-09);FilterBank.push_back(1.03470453927e-08);FilterBank.push_back(-1.43632948779e-09);FilterBank.push_back(-1.34919775398e-09);FilterBank.push_back(5.26113255735e-10);FilterBank.push_back(6.73233649018e-11);FilterBank.push_back(-8.27825652253e-11);FilterBank.push_back(1.10169293459e-11);FilterBank.push_back(6.29153731703e-12);FilterBank.push_back(-2.48478923756e-12);FilterBank.push_back(2.62649650406e-14);FilterBank.push_back(1.80866123627e-13);FilterBank.push_back(-4.24981781957e-14);FilterBank.push_back(-4.56339716212e-16);FilterBank.push_back(2.04509967678e-15);FilterBank.push_back(-4.40530704248e-16);FilterBank.push_back(4.30459683955e-17);FilterBank.push_back(-1.71615245108e-18);
                        break;
                    case sym8:
                        FilterBank.push_back(0.032223100604);FilterBank.push_back(-0.012603967262);FilterBank.push_back(-0.099219543577);FilterBank.push_back(0.297857795606);FilterBank.push_back(0.803738751807);FilterBank.push_back(0.497618667633);FilterBank.push_back(-0.029635527646);FilterBank.push_back(-0.075765714789);
                        break;
                    case sym16:
                        FilterBank.push_back(0.001889950333);FilterBank.push_back(-0.000302920515);FilterBank.push_back(-0.014952258337);FilterBank.push_back(0.003808752014);FilterBank.push_back(0.049137179674);FilterBank.push_back(-0.027219029917);FilterBank.push_back(-0.051945838108);FilterBank.push_back(0.364441894835);FilterBank.push_back(0.777185751701);FilterBank.push_back(0.481359651258);FilterBank.push_back(-0.061273359068);FilterBank.push_back(-0.143294238351);FilterBank.push_back(0.007607487325);FilterBank.push_back(0.031695087811);FilterBank.push_back(-0.000542132332);FilterBank.push_back(-0.003382415951);
                        break;
                    case coif6:
                        FilterBank.push_back(-0.072732619513);FilterBank.push_back(0.337897662458);FilterBank.push_back(0.852572020212);FilterBank.push_back(0.384864846864);FilterBank.push_back(-0.072732619513);FilterBank.push_back(-0.015655728135);
                        break;
                    case coif12:
                        FilterBank.push_back(0.016387336464);FilterBank.push_back(-0.041464936782);FilterBank.push_back(-0.067372554722);FilterBank.push_back(0.386110066823);FilterBank.push_back(0.812723635450);FilterBank.push_back(0.417005184424);FilterBank.push_back(-0.076488599079);FilterBank.push_back(-0.059434418647);FilterBank.push_back(0.023680171946);FilterBank.push_back(0.005611434819);FilterBank.push_back(-0.001823208871);FilterBank.push_back(-0.000720549445);
                        break;
                    case coif18:
                        FilterBank.push_back(-0.003793512864);FilterBank.push_back(0.007782596427);FilterBank.push_back(0.023452696142);FilterBank.push_back(-0.065771911282);FilterBank.push_back(-0.061123390003);FilterBank.push_back(0.405176902410);FilterBank.push_back(0.793777222626);FilterBank.push_back(0.428483476378);FilterBank.push_back(-0.071799821619);FilterBank.push_back(-0.082301927107);FilterBank.push_back(0.034555027573);FilterBank.push_back(0.015880544864);FilterBank.push_back(-0.009007976137);FilterBank.push_back(-0.002574517689);FilterBank.push_back(0.001117518771);FilterBank.push_back(0.000466216960);FilterBank.push_back(-0.000070983303);FilterBank.push_back(-0.000034599773);
                        break;
                    case coif24:
                        FilterBank.push_back(0.000892313669);FilterBank.push_back(-0.001629492013);FilterBank.push_back(-0.007346166328);FilterBank.push_back(0.016068943965);FilterBank.push_back(0.026682300156);FilterBank.push_back(-0.081266699681);FilterBank.push_back(-0.056077313317);FilterBank.push_back(0.415308407030);FilterBank.push_back(0.782238930921);FilterBank.push_back(0.434386056491);FilterBank.push_back(-0.066627474263);FilterBank.push_back(-0.096220442034);FilterBank.push_back(0.039334427123);FilterBank.push_back(0.025082261845);FilterBank.push_back(-0.015211731528);FilterBank.push_back(-0.005658286687);FilterBank.push_back(0.003751436157);FilterBank.push_back(0.001266561929);FilterBank.push_back(-0.000589020756);FilterBank.push_back(-0.000259974552);FilterBank.push_back(0.000062339034);FilterBank.push_back(0.000031229876);FilterBank.push_back(-0.000003259680);FilterBank.push_back(-0.000001784985);
                        break;
                    case coif30:
                        FilterBank.push_back(-0.000212080840);FilterBank.push_back(0.000358589688);FilterBank.push_back(0.002178236358);FilterBank.push_back(-0.004159358782);FilterBank.push_back(-0.010131117521);FilterBank.push_back(0.023408156788);FilterBank.push_back(0.028168028974);FilterBank.push_back(-0.091920010569);FilterBank.push_back(-0.052043163181);FilterBank.push_back(0.421566206733);FilterBank.push_back(0.774289603730);FilterBank.push_back(0.437991626216);FilterBank.push_back(-0.062035963969);FilterBank.push_back(-0.105574208714);FilterBank.push_back(0.041289208754);FilterBank.push_back(0.032683574270);FilterBank.push_back(-0.019761778945);FilterBank.push_back(-0.009164231163);FilterBank.push_back(0.006764185449);FilterBank.push_back(0.002433373213);FilterBank.push_back(-0.001662863702);FilterBank.push_back(-0.000638131343);FilterBank.push_back(0.000302259582);FilterBank.push_back(0.000140541150);FilterBank.push_back(-0.000041340432);FilterBank.push_back(-0.000021315027);FilterBank.push_back(0.000003734655);FilterBank.push_back(0.000002063762);FilterBank.push_back(-0.000000167443);FilterBank.push_back(-0.000000095177);
                        break;
                    case beylkin18:
                        FilterBank.push_back(0.099305765374353);FilterBank.push_back(0.424215360812961);FilterBank.push_back(0.699825214056600);FilterBank.push_back(0.449718251149468);FilterBank.push_back(-0.110927598348234);FilterBank.push_back(-0.264497231446384);FilterBank.push_back(0.026900308803690);FilterBank.push_back(0.155538731877093);FilterBank.push_back(-0.017520746266529);FilterBank.push_back(-0.088543630622924);FilterBank.push_back(0.019679866044322);FilterBank.push_back(0.042916387274192);FilterBank.push_back(-0.017460408696028);FilterBank.push_back(-0.014365807968852);FilterBank.push_back(0.010040411844631);FilterBank.push_back(0.0014842347824723);FilterBank.push_back(-0.002736031626258);FilterBank.push_back(0.0006404853285212);
                        break;
                    case vaidyanathan24:
                        FilterBank.push_back(-0.000062906118);FilterBank.push_back(0.000343631905);FilterBank.push_back(-0.000453956620);FilterBank.push_back(-0.000944897136);FilterBank.push_back(0.002843834547);FilterBank.push_back(0.000708137504);FilterBank.push_back(-0.008839103409);FilterBank.push_back(0.003153847056);FilterBank.push_back(0.019687215010);FilterBank.push_back(-0.014853448005);FilterBank.push_back(-0.035470398607);FilterBank.push_back(0.038742619293);FilterBank.push_back(0.055892523691);FilterBank.push_back(-0.077709750902);FilterBank.push_back(-0.083928884366);FilterBank.push_back(0.131971661417);FilterBank.push_back(0.135084227129);FilterBank.push_back(-0.194450471766);FilterBank.push_back(-0.263494802488);FilterBank.push_back(0.201612161775);FilterBank.push_back(0.635601059872);FilterBank.push_back(0.572797793211);FilterBank.push_back(0.250184129505);FilterBank.push_back(0.045799334111);
                        break;
                    default:
                        return false;
                }
            return true;
            }
void DWT::Transformada_Wavelet(double* f, long n,bool ordem){
                    int ch=FilterBank.size();

                    double *g=new double[FilterBank.size()];
                    double *h=new double[FilterBank.size()];

                    for(int i=0;i<FilterBank.size();i++)
                        h[i]=FilterBank[i];
                    for(int i=0;i<FilterBank.size();i++){
                        g[i]=h[ch-i-1];
                        if(i%2!=0)
                            g[i]*=-1;
                    }
                    int cg=ch;
                    long j=0;
                    double* t = new double[n];
                    if(ordem) // n de normal para wavelet
                        {
                        for(long i=0;i<n;i+=2)  //trend
                            {
                            t[j]=0;
                            for(long k=0;k<ch;k++)
                                t[j]+=f[(i+k)%n]*h[k];
                            j++;
                        }
                        for(long i=0;i<n;i+=2) //fluctuation
                            {
                            t[j]=0;
                            for(long k=0;k<cg;k++)
                                t[j]+=f[(i+k)%n]*g[k];
                            j++;
                        }
                    }
                    else // i de invertido para wavelet packet
                        {
                        for(long i=0;i<n;i+=2) //fluctuation
                            {
                            t[j]=0;
                            for(long k=0;k<cg;k++)
                                t[j]+=f[(i+k)%n]*g[k];
                            j++;
                        }
                        for(long i=0;i<n;i+=2)  //trend
                            {
                            t[j]=0;
                            for(long k=0;k<ch;k++)
                                t[j]+=f[(i+k)%n]*h[k];
                            j++;
                        }
                    }
                    for(long i=0;i<n;i++)
                        f[i]=t[i];

                delete(g);
                delete(h);
                delete(t);
            }

