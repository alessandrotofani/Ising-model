
#include <iostream>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <iomanip>
#include <TMath.h>
#include <ctime>
#include <time.h>
#include <cmath>

using namespace std;

const int ms = 10000; // dimensione massima del lato


int avanti(int posizione, int lato) {
	posizione = posizione + 1;
	if (posizione == lato)
		posizione = 0;
	return posizione;
}

int indietro(int posizione, int lato) {
	posizione = posizione - 1;
	if (posizione == -1)
		posizione = 0;
	return posizione;
}


double** set_spin(double **s,double lato) { //setta gli spin in modo random

  for (int i = 0; i < lato; i++) {
		for (int j = 0; j < lato; j++) {

		    s[i][j] = rand()%2;
		    if(s[i][j] == 0){ s[i][j] = -1;}
		}
	}
	return s;
}

bool** set_spin_flippato(bool **spin_flippato,double lato) {
	for (int i = 0; i < lato; i++) {
		for (int j = 0; j < lato; j++) {

			spin_flippato[i][j] = false;
		}
	}
	return spin_flippato;

}

double* binning(double* magnetizzazione, int &iterazioni){

  int iterazioni_new = iterazioni / 2;
  double* magnetizzazione_new = new double [iterazioni_new];
  for(int i = 0; i < iterazioni_new; i++){
    magnetizzazione_new[i] = (magnetizzazione[2*i] + magnetizzazione[2*i+1])/2;}
  iterazioni = iterazioni/2;
  for(int i = 0; i<iterazioni; i++){
    magnetizzazione[i] = magnetizzazione_new[i];}

    return magnetizzazione;

}

double media(double* x, double misure){  //calcola la media di x su un numero di misure dato
  double valore_medio = 0;
  for(int i = 0; i<misure;i++){
    valore_medio += x[i]/misure;
  }
  return valore_medio;}

double varianza_(double* x, double x_medio, double misure){ //calcola la varianza sulla media di x rispetto a x medio su un numero di misure dato
  double varianza_x = 0;
  double somma_quadrati = 0;

  for(int i = 0; i<misure; i++){
    somma_quadrati += pow(x[i],2);}
  for(int i = 0; i<misure; i++){
    varianza_x = (somma_quadrati/ misure - pow(x_medio,2))/misure;
  }
  return varianza_x;}


void check_magnetizzazione(double magnetizzazione_finale, double counter){
  if( magnetizzazione_finale > 1 || magnetizzazione_finale<-1){ counter = counter -1;} }

void check_varianza(double &varianza){
  if(varianza<0){
    varianza = -varianza;}}
void check_positivo(double &x){
  if(x<0){ x= -x;}
}

//ciclo  creato nell'heap
// passo alla funzione gli indirizzi di Ed e m per poterli modificare
double** Ciclo_Monte_Carlo(double **s, bool** spin_flippato, int lato, double beta, double h, double &Ed, double &m){

	double M_totale = 0;
	double E_totale = 0;
	double spin_corrente = 0;
	double somma_vicini = 0;
	double forza_efficace = 0;
	double new_spin = 0;
	double random = 0;
	int avanti_x, avanti_y, indietro_x, indietro_y;
	int siti_visitati = 0;
	int i = 0;
	int j = 0;

	while(siti_visitati<pow(lato,2)){
	  i = rand()%lato;
	  j = rand()%lato;
	  if( spin_flippato[i][j] == false){
			spin_corrente = s[i][j];
			avanti_x = avanti(i, lato);
			avanti_y = avanti(j, lato);
			indietro_x = indietro(i, lato);
			indietro_y = indietro(j, lato);
			if (avanti_x == 0) {
				s[avanti_x][j] = 0;
			}
			if (avanti_y == 0) {
				s[i][avanti_y] = 0;
			}
			if (indietro_x == 0) {
				s[indietro_x][j] = 0;
			}
			if (indietro_y == 0) {
				s[i][indietro_y] = 0;
			}
			somma_vicini = s[i][avanti_y] + s[i][indietro_y] + s[avanti_x][j] + s[indietro_x][j];
			forza_efficace = somma_vicini + h;
			random =((double) rand()/RAND_MAX); // random tra 0 e 1
			if (exp(-(beta * forza_efficace * spin_corrente * 2)) > random) {
				new_spin = -spin_corrente;
			}
			else new_spin = spin_corrente;
			s[i][j] = new_spin;
			M_totale = M_totale + new_spin;
			E_totale = E_totale -(0.5 * somma_vicini + h) * new_spin;
			m = M_totale / pow(lato, 2);
			Ed = E_totale / pow(lato, 2);
			siti_visitati++;
			spin_flippato[i][j] = true;
		}

	}






	return s;
}

void Fit()
{

  const  int pmisure = 20; // sono le magnetizzazioni finali che trova
  const int bin = 8;
  bool binned = false; // se true, fa il binning
  int iterazioni = 100000; // numero di volte che ripete il ciclo monte carlo //300
  int iterazioni_iniziali = iterazioni;
  int lato = 30; //con 30 o 40 viene
  double Iterazioni = iterazioni; //serve quando calcolo le medie
	int random = 0;
	double beta = 0.5;
	//double beta_finale[pmisure]; // serve per fare il grafico
	double sbeta_finale[pmisure];
	double h = 0;
	double Ed = 0;
	double m = 0;
	int seed = 0;
	int iterazione = 0;
	double* store_m = new double [iterazioni];
	double* magnetizzazione = new double [iterazioni];
	double* Energyd = new double [iterazioni];
	double magnetizzazione_finale[pmisure]; //media delle magnetizzazioni dopo le iterazioni//serve per fare il grafico Magnetizzazione(beta)
	double smagnetizzazione_finale[pmisure];
	double Energyd_finale[pmisure];
	double sEnergyd_finale[pmisure];
	double Energia[pmisure];
	double sEnergia[pmisure];
	double Pmisure = pmisure;
	double std[pmisure];
	double varianza[pmisure];
	double std_Energyd[pmisure];
	double varianza_Energyd[pmisure];
	int iterazioni_new = iterazioni / 2;
        double* magnetizzazione_new = new double [iterazioni_new];

	//autocorrelazione
	double t_auto[pmisure]; //tempo di autocorrelazione
	double step_1[pmisure];
	double t_auto_energia[pmisure]; //tempo di autocorrelazione
	double step_1_energia[pmisure];
	double errore_magnetizzazione[pmisure];
	double errore_energia[pmisure];
	double varianza_corretta[pmisure];
	double varianza_corretta_energia[pmisure];

	//capacita' termica
	double capacita_termica[pmisure];
	double scapacita_termica[pmisure];

	//suscettibiita' magnetiza
	double X[pmisure];
	double sX[pmisure];
	double Magnetizzazione[pmisure];
	double sMagnetizzazione[pmisure];


	srand (time(NULL)); //seed del random generator

	double beta_finale[] = {0.1, 0.2, 0.3, 0.4, 0.42, 0.43, 0.44 , 0.45 , 0.46, 0.47,  0.48, 0.49, 0.5, 0.52, 0.55,  0.6, 0.7, 0.8, 0.9, 1}; //20 misure
	  //beta_finale[] = {0.1, 0.3, 0.4, 0.44, 0.5, 0.54, 0.58, 0.6, 0.8, 1}; //10 misure


	for(int i = 0; i<pmisure; i++){
	  magnetizzazione[i] = 0;
	  Energyd[i] = 0;
	  magnetizzazione_finale[i] = 0;
	  smagnetizzazione_finale[i] = 0;
	  Energyd_finale[i] = 0;
	  sEnergyd_finale[i] = 0;
	  Energia[i] = 0;
	  sEnergia[i] = 0;
	  std[i] = 0;
	  varianza[i] = 0;
	  std_Energyd[i] = 0;
	  varianza_Energyd[i] = 0;
	  t_auto[i] = 0;
	  step_1[i] = 0;
	  t_auto_energia[i] = 0;
	  step_1_energia[i] = 0;
	  errore_magnetizzazione[i] = 0;
	  errore_energia[i] = 0;
	  varianza_corretta[i] = 0;
	  varianza_corretta_energia[i] = 0;
	  store_m[i] = 0;
	  capacita_termica[i] = 0;
	  scapacita_termica[i] = 0;
	  X[i] = 0;
	  sX[i] = 0;
	  Magnetizzazione[i] = 0;
	  sMagnetizzazione[i] = 0;
	}

	for(int i = 0; i<iterazioni; i++){  //inizializzo il vettore
	  magnetizzazione[i] = 0;
	  Energyd[i] = 0;
	}


  	double** s; //contiene una configurazione del sistema
	s = new double* [lato];
	for (int i = 0; i < lato; i++) {
		s[i] = new double[lato];
	}

	bool** spin_flippato;
	spin_flippato = new bool* [lato];
	for (int i = 0; i < lato; i++) {
		spin_flippato[i] = new bool[lato];
	}

	set_spin(s,lato);
	set_spin_flippato(spin_flippato,lato);

	for(int counter = 0; counter < pmisure; counter++){
	  beta = beta_finale[counter];
	  cout<<" beta = "<<beta<<endl;
	  iterazioni = iterazioni_iniziali;
	for (int k = 0; k < iterazioni; k++) {

	  Ciclo_Monte_Carlo(s, spin_flippato, lato, beta, h, Ed, m);
	  set_spin_flippato(spin_flippato,lato);
		store_m[k] = m;
		while( m > 1 || m < -1){
		  m = store_m[k-1];
		  Ciclo_Monte_Carlo(s, spin_flippato, lato, beta, h, Ed, m);}
		if(m>0){ magnetizzazione[k] = m;}
		    else magnetizzazione[k] = -m;
		  Energyd[k] = Ed;
		}

	if(binned){
	for( int i = 0; i<bin; i++){
	  binning(magnetizzazione,iterazioni);
	  Iterazioni = iterazioni;
	}}


	magnetizzazione_finale[counter] = media(magnetizzazione, Iterazioni);
	varianza[counter] = varianza_(magnetizzazione, magnetizzazione_finale[counter], Iterazioni); //varianza della media
	check_varianza(varianza[counter]);
	std[counter] = sqrt(varianza[counter]);


	 Energyd_finale[counter] +=media(Energyd, Iterazioni);
	 varianza_Energyd[counter] += varianza_(Energyd, Energyd_finale[counter], Iterazioni); //varianza della media
	 check_varianza(varianza_Energyd[counter]);
	 std_Energyd[counter] = sqrt(varianza_Energyd[counter]);

	 //calcolo tempo di autocorrelazione integrato

	 for( int t = 0; t<iterazioni; t++){
	    for(int i = 0; i<iterazioni;i++){
	      if(i+t<iterazioni){
		step_1[counter] += (magnetizzazione[i]*magnetizzazione[i+t] - Iterazioni*pow(magnetizzazione_finale[counter],2))/Iterazioni;
		step_1_energia[counter] += (Energyd[i]*Energyd[i+t] - Iterazioni*pow(Energyd_finale[counter],2))/Iterazioni;

	      }}}


	 t_auto[counter] = 0.5 * (step_1[counter])/(varianza[counter]); //tempo di autocorrelazione integrato
	 t_auto_energia[counter] = 0.5 * (step_1_energia[counter])/(varianza_Energyd[counter]); //tempo di autocorrelazione integrato
	 check_positivo( t_auto[counter]);
	 check_positivo( t_auto_energia[counter]);
	 //calcolo l'errore corretto sul valore medio
	 for(int t = 0; t<iterazioni; t++){
	   for( int i = 0; i<iterazioni; i++){
	      if( i+t<iterazioni){
		varianza_corretta[counter] += (1/pow(Iterazioni,2)) * (Iterazioni - abs(t))*(magnetizzazione[i]*magnetizzazione[i+t] - magnetizzazione_finale[counter])/Iterazioni;
		varianza_corretta_energia[counter] += (1/pow(Iterazioni,2)) * (Iterazioni - abs(t))*(Energyd[i]*Energyd[i+t] - Energyd_finale[counter])/Iterazioni;
	      }}}

	 check_varianza(varianza_corretta[counter]);
	 check_varianza(varianza_corretta_energia[counter]);

	 //errore corretto sui valori medi
	 errore_magnetizzazione[counter] = sqrt(varianza_corretta[counter]);
	 errore_energia[counter] = sqrt(varianza_Energyd[counter]*2*t_auto_energia[counter]/Iterazioni); //sqrt(varianza_corretta_energia[counter]);

	 //stampo i risultati
      	 cout<<" magnetizzazione = "<<magnetizzazione_finale[counter]<<endl;
	 cout<<" smagnetizzazione = "<<std[counter]<<endl;
	 cout<<" Densita' di energia = "<<Energyd_finale[counter]<<endl;
	 cout<<" sDensita' di energia = "<<std_Energyd[counter]<<endl;
	 cout<<" tempo di autocorrelazione = "<<t_auto[counter]<<endl;
	 cout<<" tempo di autocorrelazione energia = "<<t_auto_energia[counter]<<endl;
	 check_magnetizzazione(magnetizzazione_finale[counter],counter);


	} //fine del ciclo for con counter -> pmisure



  for(int i = 0; i<pmisure; i++){
    Energia[i] =( Energyd_finale[i] +2 ) * pow(lato,2) ; // metto il +2 per scalare tutte le energie e mettere lo zero corrispondente all'energia minore
    sEnergia[i] = std_Energyd[i] * pow(lato,2);
    Magnetizzazione[i] = magnetizzazione_finale[i] * pow(lato,2);
    sMagnetizzazione[i] = std[i] * pow(lato,2);

  }

  //calcolo la capacita' termica del sistema

  for(int i=0; i<pmisure;i++){
    capacita_termica[i] = pow(sEnergia[i]*beta_finale[i],2);
    scapacita_termica[i] = pow(sEnergia[i]/sqrt(2*(Iterazioni-1)) , 2);
    X[i] =  pow(sMagnetizzazione[i]*beta_finale[i],2);
    sX[i] = pow(sMagnetizzazione[i]/sqrt(2*(Iterazioni-1)) , 2);
    cout<<"Beta = "<<beta_finale[i]<<endl;
    cout<<" Capacita' termica =( "<<capacita_termica[i]<<" +- "<<scapacita_termica[i]<<" )"<<endl;
    cout<<" Suscettibilita' magnetica  =( "<<X[i]<<" +- "<<sX[i]<<" )"<<endl;

  }

	for(int i = 0; i<pmisure; i++){
	  sbeta_finale[i] = 0;}

 cout << "\n\n --- Relazione tra beta e densita' di magnetizzazione ---" <<endl;

  //  Grafico 1
  TCanvas *c1 = new TCanvas("c1","Densita' di magnetizzazione(beta)",200,10,600,400);
  //c1->SetLogx();
  c1->SetFillColor(0);
  c1->cd();
  TGraphErrors *g1 = new TGraphErrors(pmisure,beta_finale,magnetizzazione_finale,sbeta_finale,std);
  g1->SetMarkerSize(0.6);
  g1->SetMarkerStyle(21);
  g1->SetTitle("Densita' di magnetizzazione(beta)");
  g1->GetXaxis()->SetTitle("beta");
  g1->GetYaxis()->SetTitle("Densita' di magnetizzazione ");
  g1->Draw("ACP");


  // Fit 1
  cout << "\n\n Ipotesi funzione: pow([1]-x,1/8) \n" <<endl;
  TF1 *funz1 = new TF1("funz1","pow((1/[0] - 1/x), [1]) ",0.5,0.65); //0.5 funziona
  funz1->SetParameter(0,0.48);
  funz1->SetParameter(1,0.120);
  funz1->SetLineColor(1);
  g1->Fit(funz1,"RM+");
  cout<<endl;
  float beta_critico = funz1->GetParameter(0);
  float sbeta_critico = funz1->GetParError(0);
  double beta_critico_teorico = 1/2.27;
  cout<<"Beta critico =(  "<<beta_critico<<" +- "<<sbeta_critico<<" )"<<endl;
  cout<<"Beta critico teorico = "<<beta_critico_teorico<<endl;


  float Esponente_b = funz1->GetParameter(1);
  float sEsponente_b = funz1->GetParError(1);
  double esponente_b_teorico = 0.125;
  cout<<"Esponente critico =(  "<<Esponente_b<<" +- "<<sEsponente_b<<" )"<<endl;
  cout<<"Esponente critico teorico = "<<esponente_b_teorico<<endl;


  cout << "Chi^2:" << funz1->GetChisquare() << ", number of DoF: " << funz1->GetNDF() << " (Probability: " << funz1->GetProb() << ")." << endl;

  //relazione tra Energia e beta




   //  Grafico 2
  TCanvas *c2 = new TCanvas("c2","Energia(beta)",200,10,600,400);
  c2->SetFillColor(0);
  c2->cd();
  TGraphErrors *g2 = new TGraphErrors(pmisure,beta_finale,Energia,sbeta_finale,sEnergia);
  g2->SetMarkerSize(0.6);
  g2->SetMarkerStyle(21);
  g2->SetTitle("Energia(beta)");
  g2->GetXaxis()->SetTitle("beta");
  g2->GetYaxis()->SetTitle("Energia");
  g2->Draw("ACP");

  //Cpacita' Termica
  cout<<"Capacita' termica in funzione di beta"<<endl;

   //  Grafico 3
  TCanvas *c3 = new TCanvas("c3","Capacita' termica(beta)",200,10,600,400);
  c3->SetFillColor(0);
  c3->cd();
  TGraphErrors *g3 = new TGraphErrors(pmisure,beta_finale,capacita_termica,sbeta_finale,scapacita_termica);
  g3->SetMarkerSize(0.6);
  g3->SetMarkerStyle(21);
  g3->SetTitle("Capacita' termica(beta)");
  g3->GetXaxis()->SetTitle("beta");
  g3->GetYaxis()->SetTitle("Capacita' termica ");
  g3->Draw("ACP");



  //Suscettivita' magnetica
  cout<<endl;
  cout<<"Suscettivita' magnetica in funzione di beta"<<endl;

 //  Grafico 4
  TCanvas *c4 = new TCanvas("c4","Suscettivita' magnetica(beta)",200,10,600,400);
  c4->SetFillColor(0);
  c4->cd();
  TGraphErrors *g4 = new TGraphErrors(pmisure,beta_finale,X,sbeta_finale,sX);
  g4->SetMarkerSize(0.6);
  g4->SetMarkerStyle(21);
  g4->SetTitle("X(beta)");
  g4->GetXaxis()->SetTitle("beta");
  g4->GetYaxis()->SetTitle("X ");
  g4->Draw("ACP");


  // Fit 4
  cout << "\n\n Ipotesi funzione: [0]-pow(abs((-1/0.45 + 1/x)), [1]) \n" <<endl;
  TF1 *funz4 = new TF1("funz4","[0]-pow(abs((-1/[2] + 1/x)), [1]) ",0.4,0.54);
  funz4->SetParameter(0,0.8);
  funz4->SetParameter(1,1.75);
  funz4->SetParameter(2,0.46);
  funz4->SetLineColor(1);
  g4->Fit(funz4,"RM+");
  cout<<endl;

  float gamma= funz4->GetParameter(1);
  float sgamma = funz4->GetParError(1);
  double gamma_teorico = 1.75;
  cout<<"Esponente critico =(  "<<gamma<<" +- "<<sgamma<<" )"<<endl;
  cout<<"Esponente critico teorico = "<<gamma_teorico<<endl;


  cout << "Chi^2:" << funz4->GetChisquare() << ", number of DoF: " << funz4->GetNDF() << " (Probability: " << funz4->GetProb() << ")." << endl;


  }
