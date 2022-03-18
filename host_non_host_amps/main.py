def main():
    import os
    from utils.download_files import inputsgen
    from utils.load_data import input_info
    from utils.distributions import distributions
    from utils.ampsets_generator import ampsets
    from utils.boxstats import plot_test
    
    os.makedirs('data/', exist_ok=True)
    os.makedirs('analysis/', exist_ok=True)

    print('Download inputs')
    inputsgen()
    print('Load data')
    data, spheres = input_info()
    print('Calculate distributions')
    distributions(data, spheres)
    print('Select sets and normalize')
    df = ampsets(data)
    print('Plot AMPs per environment corrected')
    plot_test(df)


if __name__ == '__main__':
    main()
        
