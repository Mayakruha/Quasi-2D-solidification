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
    Q,DeltaR,SpUs:array of Extended;{Basic var}
    nk:array of word;{Basic var}
    ksiC:array[0..2] of Extended;
    m,count,n0:word;{Basic var}
    F:TextFile;
    Tsol,Tlik,T0,beta,Hl,KsiI:Extended;{Key parametrs}
    dr,dtau,dZ,Value,Tmin,KsiZ:Extended;{Basic var}
    {Local var}
    Hk,Tk,dH,ValueLoc,eps:Extended;
    dx,Nf0,Nf1,Nf2:array[0..1] of Extended;
    Flag:boolean;
  public
    { Public declarations }
    function HeatFlow(z,Ts:Extended):Extended;
    function Cef(Temp:Extended):Extended;
    function FuncTemp(Value:Extended):Extended;
    function Temperature(Value:Extended):Extended;
    procedure FindKsiC(var ksi:array of Extended;SpUs0,KsiZ0:Extended;j0,i0:word);
    procedure NormalForce(var Nf:array of Extended;SpUs0,KsiZ0:Extended; j0:word);
    procedure SpeedUsadka(var SpUs0,KsiZ0:Extended;j0:word);
  end;

const
{----------chemical compound of steel--------------------------}
  C_C=0.12;
  C_Mn=1.3;
  C_Si=0.8;
  C_P=0.035;
  C_S=0.04;
  C_Al=0.0;
  C_Cu=0.3;
  C_Ni=0.3;
  C_Cr=0.3;
{-----------steel parametr----------------------}
  lamda=29{Vt/m*K}; {теплопроводность слитка}
  L=272{kJ/kg}*1000;{теплота кристаллизации}
  ro=7200{kg/m3};{плотность стали}
  Cl=500{J/kg*K};{теплоемкость жидкой стали}
  Cr=680{J/kg*K};{теплоемкость твердой стали}
  k=0.1454;{степень в законе ползучести стали}
  Tc=265.96;{K}{коэффициент в законе ползучести стали}
  Ac=11348;{MPa}{коэффициент в законе ползучести стали}

{---------- casting parametr---------------------------}
  Lz=1200{mm}/1000;{расчетная длина}
  R=314{mm}/1000;{радиус слитка}
  v=0.28{m/min}; {скорость разливки}
  dTemp=40;{температура перегрева}
{---------- cooling parametr-----------------------------}
  Tsr=30;{температура окружающей среды, Цельсий}
{-------- calc parametr-------------------------}
  n=300;
  kj=1; {коэффициент сходимости температурной задачи}
  Epsilon=0.0001;

var
  Form1: TForm1;

implementation

Uses Math, Matchad;
{$R *.DFM}
{--------------- процедуры для решения температурной задачи ------------}
function TForm1.HeatFlow(z,Ts:Extended):Extended;
begin
  if z<0.64 then Result:=(1.236*exp(-1.6*z/v)+0.364)*1000000 {кристаллизатор}
 { else if z<0.96 then Result:=4.5*(Power(Ts/100,4)-Power(Tsr/100,4))+40*(Ts-Tsr){1-ая зона, 1-ый участок}
  else if z<1.3 then Result:=500*(Ts-Tsr){1-ая зона, 2-ой участок}
  {else if z<5.5 then Result:=800*(Ts-Tsr){2-ая и 3-я зоны}
  {else Result:=4.5*(Power(Ts/100,4)-Power(Tsr/100,4))+40*(Ts-Tsr)};{рольганг}
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

{--------------- end: процедуры для решения температурной задачи ------------}

{--------------- процедуры для расчета усадки ------------}

procedure TForm1.FindKsiC(var ksi:array of Extended;SpUs0,KsiZ0:Extended;j0,i0:word);
var i:word;
begin
  ksi[0]:=2*beta*SpT[j0*(n+1)+i0]-KsiZ0/2-(KsiZ0*R*R/2+SpUs0*R)/i0/i0/dr/dr;
  ksi[1]:=(KsiZ0*R*R/2+SpUs0*R)/i0/i0/dr/dr-beta*SpT[j0*(n+1)+i0]-KsiZ0/2;
  ksi[2]:=KsiZ0-beta*SpT[j0*(n+1)+i0];
  for i:=i0 to n-1 do
    begin
      ValueLoc:=3*beta*(SpT[j0*(n+1)+i]+SpT[j0*(n+1)+i+1])*(i+0.5)/i0/i0/2;
      ksi[0]:=ksi[0]+ValueLoc;
      ksi[1]:=ksi[1]-ValueLoc;
    end;
end;

procedure TForm1.NormalForce(var Nf:array of Extended;SpUs0,KsiZ0:Extended;j0:word);
var
  i:word;
begin
  Nf[0]:=0;
  Nf[1]:=0;
  for i:=nk[j0] to n do
    begin
      FindKsiC(ksiC,SpUs0,KsiZ0,j0,i);
      KsiI:=Sqrt(2*((ksiC[0]-ksiC[1])*(ksiC[0]-ksiC[1])+
          (ksiC[0]-ksiC[2])*(ksiC[0]-ksiC[2])+
          (ksiC[1]-ksiC[2])*(ksiC[1]-ksiC[2])))/3;
      if KsiI=0 then ValueLoc:=0
      else ValueLoc:=2*Ac*Power(KsiI,k-1)*exp(-(T[j0*(n+1)+i]+273)/Tc)/3;
      if (i=nk[j0])or(i=n) then
        begin
          Nf[0]:=Nf[0]+ValueLoc*(ksiC[1]-ksiC[0])*dr/2;
          Nf[1]:=Nf[1]+ValueLoc*(ksiC[2]-ksiC[0])*dr/2;
        end
      else
        begin
          Nf[0]:=Nf[0]+ValueLoc*(ksiC[1]-ksiC[0])*dr;
          Nf[1]:=Nf[1]+ValueLoc*(ksiC[2]-ksiC[0])*dr;
        end;
    end;
end;

procedure TForm1.SpeedUsadka(var SpUs0,KsiZ0:Extended;j0:word);
begin
  if SpUs0=0 then dx[0]:=-beta*100*R
  else dx[0]:=-100*Epsilon*SpUs0;
  if KsiZ0=0 then dx[1]:=-beta*100
  else dx[1]:=-100*Epsilon*KsiZ0;
  repeat
    NormalForce(Nf0,SpUs0,KsiZ0,j0);
    NormalForce(Nf1,SpUs0+dx[0],KsiZ0,j0);
    NormalForce(Nf2,SpUs0,KsiZ0+dx[1],j0);

    if (Nf0[0]*Nf1[0])<0 then dx[0]:=dx[0]/2
    else if ((Nf1[0]-Nf0[0])*Nf0[0])>0 then dx[0]:=-dx[0];
    if (Nf0[1]*Nf2[1])<0 then dx[1]:=dx[1]/2
    else if ((Nf2[1]-Nf0[1])*Nf0[1])>0 then dx[1]:=-dx[1];

    eps:=abs(2*dx[0]/(abs(SpUs0)+abs(dx[0])));

    if (Nf0[0]*Nf1[0])<0 then SpUs0:=SpUs0+dx[0]/2
    else if ((Nf1[0]-Nf0[0])*Nf0[0])>=0 then SpUs0:=SpUs0
    else SpUs0:=SpUs0+dx[0];
    if (Nf0[1]*Nf2[1])<0 then KsiZ0:=KsiZ0+dx[1]/2
    else if ((Nf2[1]-Nf0[1])*Nf0[1])>=0 then KsiZ0:=KsiZ0
    else KsiZ0:=KsiZ0+dx[1];

  until eps<Epsilon;
end;

procedure TForm1.ButtonClick(Sender: TObject);
var
  i,j:word;
begin
  Tsol:=1536-(200*C_C+16*C_Si+6*C_Mn+1.7*C_Cr+3.9*C_Ni+93*C_P+1100*C_S);{C}
  Tlik:=1536-(78*C_C+7.5*C_Si+4.9*C_Mn+1.3*C_Cr+3.1*C_Ni+4.7*C_Cu+3.6*C_Al+
        34.4*C_P+38*C_S);{C}
  T0:=Tlik+dTemp;
  beta:=(2.7-0.16*C_C+0.039*C_Mn-0.1*C_Si-0.019*C_Cr-0.016*C_Ni-0.5*C_P-0.25*C_S)/140000;{1/C}
  Hl:=ro*((Cl+Cr)*(Tlik-Tsol)/2+L){J/m3};
  {--------------------------------------}
  dr:=R/n;
  dtau:=kj*dr*dr*ro*Min(Cl,Cr)/4/lamda{sek};
  m:=2*Round(Lz/v/dtau*60);
  dZ:=Lz/m;
  dtau:=dZ/V*60;
  SetLength(T,(n+1)*(m+1));
  SetLength(H,(n+1)*(m+1));
  SetLength(SpT,(n+1)*m);
  SetLength(Q,m);
  SetLength(SpUs,m);
  SetLength(nk,m+1);
  SetLength(deltaR,m+1);
  Series1.Clear;
  Series2.Clear;
  for i:=0 to n do
    begin
      T[i]:=T0;
      H[i]:=FuncTemp(T0);
    end;
  Tmin:=T0;
  nk[0]:=n;
  SpUs[0]:=0;
  DeltaR[0]:=0;
  KsiZ:=0;
  count:=0;
  Flag:=True;
  for j:=0 to m-1 do
    begin
     {------------------блок расчета температуы--------------------------------}
      Q[j]:=HeatFlow(j*dz,T[n+(n+1)*j]);
      nk[j+1]:=n;
      for i:=0 to n do
        begin
          if i=0 then Value:=4*lamda*dtau/dr/dr*(T[i+1+(n+1)*j]-T[i+(n+1)*j])+
           H[i+(n+1)*j]
          else if i=n then Value:=2*lamda*dtau*(1-1/2/i)/dr/dr*(T[i-1+(n+1)*j]-T[i+(n+1)*j])
           +H[i+(n+1)*j]-2*dtau/dr*Q[j]
          else Value:=lamda*dtau/dr/dr*((1-1/2/i)*(T[i-1+(n+1)*j]-T[i+(n+1)*j])+
           (1+1/2/i)*(T[i+1+(n+1)*j]-T[i+(n+1)*j]))+H[i+(n+1)*j];

          H[i+(n+1)*(j+1)]:=Value;
          Value:=Temperature(Value);
          T[i+(n+1)*(j+1)]:=Value;
          if (nk[j+1]=n)and(Value<Tsol)then nk[j+1]:=i;
          SpT[i+(n+1)*j]:=(Value-T[i+(n+1)*j])/dtau;
          if Tmin>Value then Tmin:=Value;
        end;
      {----------------блок расчет усадки-----------------------------}
      if T[n-1+(n+1)*j]<Tsol then
        begin
          SpUs[j]:=SpUs[j-1];
          SpeedUsadka(SpUs[j],KsiZ,j);
        end
      else SpUs[j]:=0;
      deltaR[j+1]:=deltaR[j]+SpUs[j]*dtau*1000;
      {---------------------------вывод данных--------------------------------------}
      if Flag and (T[(n+1)*j]<Tsol) then
        begin
          Flag:=False;
          Label1.Caption:='Длина лунки, мм: '+ FloatToStr(j*dz);
        end;
      if count=1 then
        begin
          Series1.AddXY(j*dz*1000,deltaR[j+1]);
          Series2.AddXY(j*dz*1000,T[n+(n+1)*j]);
          {Series3.AddXY(j*dz*1000,lamda*(T[n-1+(n+1)*j]-T[n+(n+1)*j])/dr);}
          ProgressBar.Position:=Round(ProgressBar.Max*j/m);
          Application.ProcessMessages;
          count:=0;
        end;
        count:=count+1;
    end;
  {--------- Out Result--------------------}
 { AssignFile(F,'Temp.mke');
  Rewrite(F);
  Writeln(F,'температура по слябу');
  Writeln(F,'Число элементов:'+IntToStr(n*m));
  Writeln(F,'Число узлов в элементе:4');
  Writeln(F,'Минимальная координата по оси Х:0');
  Writeln(F,'Максимальная координата по оси Х:'+FloatToStr(Lz));
  Writeln(F,'Минимальноя координата по оси У:'+FloatToStr(0));
  Writeln(F,'Максимальная координата по оси У:'+FloatToStr(R));
  WriTeln(F,'Минимальное значение:'+FloatToStr(Tmin));
  WriTeln(F,'Максимальное значение:'+FloatToStr(T0));
  for j:=0 to m-1 do
    for i:=0 to n-1 do
     begin
       Write(F,FloatToStr(dZ*j)+'|'+FloatToStr(dr*i)+'|'+FloatToStr(T[i+(n+1)*j])+'|');
       Write(F,FloatToStr(dZ*(j+1))+'|'+FloatToStr(dr*i)+'|'+FloatToStr(T[i+(n+1)*(j+1)])+'|');
       Write(F,FloatToStr(dZ*(j+1))+'|'+FloatToStr(dr*(i+1))+'|'+FloatToStr(T[i+1+(n+1)*(j+1)])+'|');
       Writeln(F,FloatToStr(dZ*j)+'|'+FloatToStr(dr*(i+1))+'|'+FloatToStr(T[i+1+(n+1)*j])+'|');
     end;
  CloseFile(F);}

  AssignFile(F,'Result.txt');
  Rewrite(F);
  Writeln(F,'Расстояние от мениска, мм; Усадка, мм; Толщина корки, мм');
  for j:=0 to Round(m/10)-1 do
    begin
      Writeln(F,FloatToStr(j*10*dZ*1000)+';'+FloatToStr(-2*deltaR[j*10])+';'+FloatToStr(dr*1000*(n-nk[j*10])));
    end;
  CloseFile(F);

  {AssignFile(F,'Sigma.txt');
  Rewrite(F);
  for i:=nk[m-1] to n do
    begin
      ksi1:=ksi[m-1]-beta*SpT[(m-1)*(n+1)+i];
      if ksi1<0 then ValueLoc:=-Ac*Power(-ksi1,k)*exp(-T[(m-1)*(n+1)+i]/Tc)
      else ValueLoc:=Ac*Power(ksi1,k)*exp(-T[(m-1)*(n+1)+i]/Tc);
      Writeln(F,ToFormatMatchad((n-i)*dX)+','+ToFormatMatchad(ValueLoc));
    end;
  CloseFile(F);}

  T:=Nil;
  H:=Nil;
  SpT:=Nil;
  Q:=Nil;
  nk:=Nil;
  DeltaR:=Nil;
  ProgressBar.Position:=0;
end;

end.
