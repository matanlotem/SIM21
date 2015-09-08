"C:\Program Files (x86)\WinSCP\WinSCP.com" /command ^
    "open sftp://matanlotem@power.tau.ac.il" ^
    "option batch continue" ^
    "synchronize -delete -filemask=""|.git/; ~*"" remote C:\Users\Matan\Work\SIM21 /a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21" ^
    "keepuptodate -delete -filemask=""|.git/; ~*"" C:\Users\Matan\Work\SIM21 /a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21" ^
    exit