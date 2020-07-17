function Pause ($Message="Press any key to continue...")
{
Write-Host -NoNewLine $Message
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
Write-Host ""
}
while($true) {
cls
g++ -static -O3 -o p_z3 "prismatoid_z3/main.cpp" "prismatoid_z3/Prismatoid_Z3.cpp" -Iprismatoid_z3 -IOGDF -Lz3/lib -LOGDF/_release -lz3 -lgmp -lOGDF -lpthread 
g++ "test.cpp" -IOGDF -LOGDF/_release -lOGDF -pthread 
./a.exe
Pause
}
