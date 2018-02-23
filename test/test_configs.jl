
test_configs = [
    #
   ( "si, cubic, cluster, short",
     set_pbc!(bulk("Si", cubic=true) * 3, false),
     1.1 * rnn("Si") ),
    #
   ( "si, cubic, cluster, med",
     set_pbc!(bulk("Si", cubic=true) * 3, false),
     2.1 * rnn("Si") ),
    #
   ( "si, non-cubic cell, cluster, med",
     set_pbc!(bulk("Si") * 5, false),
     2.1 * rnn("Si") ),
    #
   ( "si, non-cubic cell, pbc",
     set_pbc!(bulk("Si") * 5, false),
     2.1 * rnn("Si") ),
    #
   ( "si, non-cubic cell, mixed bc",
     set_pbc!(bulk("Si") * 5, false),
     2.1 * rnn("Si") ),
   ]
