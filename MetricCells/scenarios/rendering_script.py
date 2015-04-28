def read_cspace_file(fname):
    f = open(fname, "r")
    fout = open( "out.txt", "w" )
    for line in f:
        info = line.split();
        x = float(info[0])
        y = float(info[1])
        z = float(info[2])
        r = float(info[3])
        w = float(info[4])
        if( w > 400 ):
            fout.write( "create_cell({0}, {1}, {2}, {3}, {4});\n".format( x, y, z, r, w ) );

    fout.flush();
    fout.close();
    f.close();

#create_cell(-2,1,0, 1, 0.5)

read_cspace_file('/Users/Yinan/Workspace/MetricCells/MetricCells/scenarios/3rarm_records.txt')
