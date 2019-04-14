unit Usadka;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls, ComCtrls, TeEngine, Series, ExtCtrls, TeeProcs, Chart;

type
  TForm1 = class(TForm)
    Button: TButton;
    ProgressBar: TProgressBar;
    Label1: TLabel;
    Chart1: TChart;
    Series1: TPointSeries;
    Chart2: TChart;
    Series2: TPointSeries;
    Series3: TPointSeries;
    procedure ButtonClick(Sender: TObject);

  private
    { Private declarations }
    T,H,SpT:array of Extended;{Basic var}
    Q,ksi,deltaB,KapaT:array of Extended;{Basic var}
    nk:array of word;{Basic var}
    F:TextFile;
    Tsol,Tlik,T0,beta,Hl,TolKorka,SigmaX:Extended;{Key parametrs}
    dX,dtau,dZ,Value,Tmin:Extended;{Basic var}
    Hk,Tk,dH,eps,ValueLoc,F0,F1,ksi1,dksi,Sum:Extended;{Local var}
    m,count:word;{Basic var}
    Flag:boolean;
  public
    { Public declarations }
    function HeatFlow(z,Ts:Extended):Extended;
    function Cef(Temp:Extended):Extended;
    function FuncTemp(Value:Extended):Extended;
    function Temperature(Value:Extended):Extended;
    function NormalForce(j,n0:word;ksi0:Extended):Extended;
    function SpeedUsadka(j,n0:word;ksi0:Extended):Extended;
  end;

const
{----------chemical compound of steel--------------------------}
  C_C=0.09;
  C_Mn=1.44;
  C_Si=0.61;
  C_P=0.012;
  C_S=0.006;
  C_Al=0.027;
  C_Cu=0.09;
  C_Ni=0.05;
  C_Cr=0.06;
{-----------steel parametr----------------------}
  lamda=29{Vt/m*K}; {теплопроводность слитка}
  L=272{kJ/kg}*1000;{теплота кристаллизации}
  ro=7200{kg/m3};{плотность стали}
  Cl=500{J/kg*K};{теплоемкость жидкой стали}
  Cr=680{J/kg*K};{теплоемкость твердой стали}
  k=0.1454;{степень в законе ползучести стали}
  Tc=266;{K}{коэффициент в законе ползучести стали}
  Ac=11348;{MPa}{коэффициент в законе ползучести стали}

{---------- casting parametr---------------------------}
  Lz=1450{mm}/1000;{расчетная длина}
  a=600{mm}/1000;{половина ширины слитка}
  b=95{mm}/1000;{половина толщины слитка}
  v=1.25{m/min}; {скорость разливки}
  dTemp=17;{температура перегрева}
{---------- cooling parametr-----------------------------}
  Tsr=300;{температура окружающей среды, Кельвин}
{-------- calc parametr-------------------------}
  n=300;
  kj=1; {коэффициент сходимости температурной задачи}
  Epsilon=0.0001;

var
  Form1: TForm1;

implementation

Uses Math, Matchad;
{$R *.DFM}

function TForm1.HeatFlow(z,Ts:Extended):Extended;
begin
  if z<0.78 then Result:=(1.155*exp(-10*z/v)+0.845)*1000000 {кристаллизатор}
 { else if z<0.96 then Result:=4.5*(Power((Ts+273)/100,4)-Power((Tsr+273)/100,4))+40*(Ts-Tsr){1-ая зона, 1-ый участок}
{  else if z<1.1 then Result:=550*(Ts-Tsr){1-ая зона, 2-ой участок}
  else if z<5.5 then Result:=1000*(Ts-Tsr){2-ая и 3-я зоны}
  else Result:=4.5*(Power((Ts+273)/100,4)-Power((Tsr+273)/100,4))+40*(Ts-Tsr);{рольганг}
end;

function TForm1.Cef(Temp:Extended):Extended;
begin
  if Temp<Tsol then Result:=Cr
  else if Temp>Tlik then Result:=Cl
  else Result:=(L+Cr*(Tlik-Temp)+Cl*(Temp-Tsol))/(Tlik-Tsol);
end;

function TForm1.FuncTemp(Value:Extended):Extended;
begin
  if Value<Tsol then Result:=(Value-Tsol)*ro*Cr
  else if Value>Tlik then Result:=(Value-Tlik)*ro*Cl+Hl
  else Result:=ro*(Value-Tsol)/(Tlik-Tsol)*(L+(Tlik-Tsol)*Cr+
      (Cl-Cr)*(Value-Tsol)/2);
end;

function TForm1.Temperature(Value:Extended):Extended;
begin
  if Value<0 then Result:=Value/ro/Cr+Tsol
  else if Value>Hl then Result:=(Value-Hl)/ro/Cl+Tlik
  else
    begin
      Tk:=Value/Hl*(Tlik-Tsol)+Tsol;
      repeat
        Hk:=FuncTemp(Tk);
        dH:=Value-Hk;
        eps:=abs(dH/Hl);
        Tk:=dH/ro/Cef(Tk)+Tk;
      until eps<Epsilon;
      Result:=Tk;
    end;
end;

function TForm1.NormalForce(j,n0:word;ksi0:Extended):Extended;
var
  i:word;
begin
  Sum:=0;
  for i:=n0 to n do
    begin
      ValueLoc:=ksi0-beta*SpT[j*(n+1)+i];
      if ValueLoc<0 then ValueLoc:=-Power(-ValueLoc,k)*exp(-(T[j*(n+1)+i]+273)/Tc)
      else ValueLoc:=Power(ValueLoc,k)*exp(-(T[j*(n+1)+i]+273)/Tc);
      if (i=n0)or(i=n) then Sum:=Sum+ValueLoc*dX/2
      else Sum:=Sum+ValueLoc*dX;
    end;
  Result:=Sum;
end;

function TForm1.SpeedUsadka(j,n0:word;ksi0:Extended):Extended;
begin
  if ksi0=0 then dksi:=-beta*100
  else dksi:=-100*Epsilon*ksi0;
  repeat
    ksi1:=ksi0+dksi;
    F0:=NormalForce(j,n0,ksi0);
    F1:=NormalForce(j,n0,ksi1);
    if (F0*F1)<0 then dksi:=dksi/2
    else if ((F1-F0)*F0)>0 then dksi:=-dksi
    else if (F0*F1)=0 then dksi:=0;
    eps:=abs(2*dksi/(abs(ksi0)+abs(ksi1)));
    if (F0*F1)<0 then ksi0:=(ksi1+ksi0)/2
    else if ((F1-F0)*F0)>=0 then ksi0:=ksi0
    else ksi0:=ksi1;
  until eps<Epsilon;
  Result:=ksi0;
end;

procedure TForm1.ButtonClick(Sender: TObject);
var
  i,j:word;
begin
  AssignFile(F,'Result.csv');
  Rewrite(F);
  Tsol:=1536-(200*C_C+16*C_Si+6*C_Mn+1.7*C_Cr+3.9*C_Ni+93*C_P+1100*C_S);{K}
  Tlik:=1536-(78*C_C+7.5*C_Si+4.9*C_Mn+1.3*C_Cr+3.1*C_Ni+4.7*C_Cu+3.6*C_Al+
        34.4*C_P+38*C_S);{K}
  T0:=Tlik+dTemp;
  beta:=(2.7-0.16*C_C+0.039*C_Mn-0.1*C_Si-0.019*C_Cr-0.016*C_Ni-0.5*C_P-0.25*C_S)/140000;{1/K}
  Hl:=ro*((Cl+Cr)*(Tlik-Tsol)/2+L){J/m3};
  {--------------------------------------}
  dX:=b/n;
  dtau:=kj*dX*dX*ro*Min(Cl,Cr)/4/lamda{sek};
  m:=2*Round(Lz/v/dtau*60);
  dZ:=Lz/m;
  dtau:=dZ/V*60;
  SetLength(T,(n+1)*(m+1));
  SetLength(H,(n+1)*(m+1));
  SetLength(SpT,(n+1)*m);
  SetLength(Q,m);
  SetLength(nk,m+1);
  SetLength(ksi,m);
  SetLength(deltaB,m+1);
  SetLength(KapaT,m+1);
  Series1.Clear;
  Series2.Clear;
  for i:=0 to n do
    begin
      T[i]:=T0;
      H[i]:=FuncTemp(T0);
    end;
  Tmin:=T0;
  SigmaX:=0;
  nk[0]:=n;
  deltaB[0]:=0;
  KapaT[0]:=0;
  count:=0;
  Flag:=True;
  for j:=0 to m-1 do
    begin
      Q[j]:=HeatFlow(j*dz,T[n+(n+1)*j]);
      nk[j+1]:=n;
      for i:=0 to n do
        begin
        {--------------------блок расчета температуры---------------------}
          if i=0 then Value:=2*lamda*dtau/dX/dX*(T[i+1+(n+1)*j]-T[i+(n+1)*j])+
           H[i+(n+1)*j]
          else if i=n then Value:=2*lamda*dtau/dX/dX*(T[i-1+(n+1)*j]-T[i+(n+1)*j])
           +H[i+(n+1)*j]-2*dtau/dX*Q[j]
          else Value:=lamda*dtau/dX/dX*(T[i-1+(n+1)*j]+
           T[i+1+(n+1)*j]-2*T[i+(n+1)*j])+H[i+(n+1)*j];
          H[i+(n+1)*(j+1)]:=Value;
          Value:=Temperature(Value);
          T[i+(n+1)*(j+1)]:=Value;
          if (nk[j+1]=n)and(Value<Tsol)then nk[j+1]:=i;
          SpT[i+(n+1)*j]:=(Value-T[i+(n+1)*j])/dtau;
          if Tmin>Value then Tmin:=Value;
        end;
      {----------------блок расчет усадки-----------------------------}
      if T[n+(n+1)*j]<Tsol then
        begin
          ksi[j]:=SpeedUsadka(j,nk[j],beta*SpT[n+(n+1)*j]);
          SigmaX:=Ac*Power(abs(2*(beta*SpT[j*(n+1)+n]-ksi[j])),k)*exp(-(T[j*(n+1)+n]+273)/Tc);
        end
      else ksi[j]:=0;
      deltaB[j+1]:=deltaB[j]+ksi[j]*a*dtau*1000;
      {---------------блок расчет скорости кривизны и толщины корки---------------------}
      KapaT[j+1]:=0;
      if (nk[j]<(n-2))and(nk[j]>2) then
        begin
          TolKorka:=((n-nk[j])*dx+(Tsol-T[nk[j]+(n+1)*j])*dx/(T[nk[j]-1+(n+1)*j]-T[nk[j]+(n+1)*j]));
          for i:=nk[j] to n do KapaT[j+1]:=KapaT[j+1]+beta*SpT[i+(n+1)*j]*(i*dX-b+TolKorka/2);
          KapaT[j+1]:=KapaT[j]+12*KapaT[j+1]*dX/TolKorka/TolKorka/Tolkorka*dtau;
        end;

      {---------------------------вывод данных--------------------------------------}
      if Flag and (T[(n+1)*j]<Tsol) then
        begin
          Flag:=False;
          Label1.Caption:='Длина лунки, мм: '+ FloatToStr(j*dz);
        end;
      if count=20 then
        begin
          Series1.AddXY(j*dz,T[j*(n+1)+n]);
          Series2.AddXY(j*dz,SigmaX);
          {Series3.AddXY(j*dz,T[(n+1)*j]);}
          Writeln(F,FloatToStr(j*dz)+';'+FloatToStr(-deltaB[j])+';'+FloatToStr(TolKorka));
          ProgressBar.Position:=Round(ProgressBar.Max*j/m);
          Application.ProcessMessages;
          count:=0;
        end;
        count:=count+1;
    end;
  {--------- Out Result--------------------}
 { AssignFile(F,'Temp.txt');
  Rewrite(F);
  Writeln(F,'температура по слябу');
  Writeln(F,'Число элементов:'+IntToStr(n*m));
  Writeln(F,'Число узлов в элементе:4');
  Writeln(F,'Минимальная координата по оси Х:0');
  Writeln(F,'Максимальная координата по оси Х:'+FloatToStr(Lz));
  Writeln(F,'Минимальноя координата по оси У:'+FloatToStr(0));
  Writeln(F,'Максимальная координата по оси У:'+FloatToStr(b));
  WriTeln(F,'Минимальное значение:'+FloatToStr(Tmin));
  WriTeln(F,'Максимальное значение:'+FloatToStr(T0));
  for j:=0 to m-1 do
    for i:=0 to n-1 do
     begin
       Write(F,FloatToStr(dZ*j)+'|'+FloatToStr(dX*i)+'|'+FloatToStr(T[i+(n+1)*j])+'|');
       Write(F,FloatToStr(dZ*(j+1))+'|'+FloatToStr(dX*i)+'|'+FloatToStr(T[i+(n+1)*(j+1)])+'|');
       Write(F,FloatToStr(dZ*(j+1))+'|'+FloatToStr(dX*(i+1))+'|'+FloatToStr(T[i+1+(n+1)*(j+1)])+'|');
       Writeln(F,FloatToStr(dZ*j)+'|'+FloatToStr(dX*(i+1))+'|'+FloatToStr(T[i+1+(n+1)*j])+'|');
     end;
  CloseFile(F); }

  CloseFile(F);

  T:=Nil;
  H:=Nil;
  SpT:=Nil;
  Q:=Nil;
  nk:=Nil;
  ksi:=Nil;
  deltaB:=Nil;
  KapaT:=Nil;
  ProgressBar.Position:=0;
end;

end.
