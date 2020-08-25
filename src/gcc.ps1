function Pause ($Message="Press any key to continue...")
{
Write-Host -NoNewLine $Message
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
Write-Host ""
}
while($true) {
cls
g++ -static -O3 "prismatoid/main.cpp" "prismatoid/prismatoid.cpp" "prismatoid/binomial.cpp" -Iprismatoid -Lz3/lib -lz3 -lgmp -lpthread -o p
g++ "test.cpp" -IOGDF -LOGDF/_release -lOGDF -pthread 
./a.exe
Pause
}
