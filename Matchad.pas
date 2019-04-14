unit Matchad;

interface

Function ToFormatMatchad(x:real):string;{}
Function IntFromText(Text:string):word;{Из текста выносит цифры}
Procedure StrokaValue(Text,Symbol:string;var x:array of string);

implementation

Uses SysUtils;

Function ToFormatMatchad(x:real):string;
var
  prom:string;
  k:byte;
begin
  prom:=FormatFloat('0.00000000000000E+',x);
  k:=Pos(',',prom);
  if k>0 then prom[k]:='.';
  Result:=prom;
end;

Function IntFromText(Text:string):word;
var
  ValueText:string;
  i:word;
begin
  ValueText:='';
  for i:=1 to Length(Text) do
    if (Text[i]>='0')and(Text[i]<='9') then ValueText:=ValueText+Text[i];
  Result:=StrToInt(ValueText);
end;

Procedure StrokaValue(Text,Symbol:string;var X:array of string);
var
  i,NValue:word;
  ValueText:string;
begin
  NValue:=0;
  ValueText:='';
  for i:=1 to Length(Text) do
    if Text[i]=Symbol then
      begin
        X[NValue]:=ValueText;
        NValue:=NValue+1;
        ValueText:='';
      end
    else if i=Length(Text) then
      begin
        if (Text[i]>='0')and(Text[i]<='9')then ValueText:=ValueText+Text[i];
        X[NValue]:=ValueText;
      end
    else if Text[i]<>' ' then  ValueText:=ValueText+Text[i];
end;

end.
