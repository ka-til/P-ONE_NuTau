def find_emin_emax(X):
    hdf = pd.HDFStore(x, mode = r)
    df = hdf.get('/I3MCTree')
    pType = df.type
    energy  = df.energy[(pType == 16 ) | (pType == -16)]
    emin = energy.min()
    emax = energy.max()
    print(emin, emax)