import argparse
from ROOT import TH1F, TCanvas, TFile, gStyle

def read_file_to_arrays(filename):
    '''
        Read the file and extract the values and errors as arrays.
    '''
    values, errors = [], []
    start_reading = False

    with open(filename, 'r') as file:
        for line in file:
            if 'dx' in line:
                dx = float(line.split()[2])
                continue
            elif 'ndata' in line:
                nx = int(line.split()[2])
                continue
            elif 'xoffset' in line:
                xoffset = float(line.split()[2])
                continue

            # Start reading after encountering '{'
            if '{' in line:
                start_reading = True
                continue

            if not start_reading:
                continue

            # Skip empty lines or lines that start with a comment character
            if not line.strip() or line.startswith('#') or line.startswith('}'):
                continue

            columns = line.split()
            values.append(float(columns[0]))
            errors.append(float(columns[1]))

    hist_specs = [nx, dx, xoffset]
    return values, errors, hist_specs

def fill_hist_from_arrays(values, errors, hist_specs, outfilepdf=None):
    '''
        Fill the histogram with the values and errors.
    '''
    nx, dx, xoffset, name, title = hist_specs
    hist = TH1F(name, title, nx, xoffset, nx * dx + xoffset)

    for ibin, (value, error) in enumerate(zip(values, errors), start=1):
        hist.SetBinContent(ibin, value)
        hist.SetBinError(ibin, error)

    if outfilepdf:
        canvas = TCanvas('canvas', 'canvas', 800, 600)
        hist.Draw('E1')
        canvas.SaveAs(outfilepdf)

    return hist

if __name__ == '__main__':

    gStyle.SetOptStat(0)
    parser = argparse.ArgumentParser(description='Plot the correlation function and same event histograms.')
    parser.add_argument('--radius', help='Radius of the particle source in fm', type=float, default=4.2)
    args = parser.parse_args()
    
    infile_CF = 'output/phemto_output_CF.dat'
    values, errors, hist_spec = read_file_to_arrays(infile_CF)
    hist_spec.extend([f'CF_{args.radius}fm', f'R = {args.radius} fm; #it{{k}}* (MeV/#it{{c}}); C(#it{{k}}*)'])
    hist_CF = fill_hist_from_arrays(values, errors, hist_spec, 'output/correlation.pdf')
    
    infile_SE = 'output/phemto_output_kstarNum.dat'
    values, errors, hist_spec = read_file_to_arrays(infile_SE)
    hist_spec.extend([f'SE_{args.radius}fm', f'R = {args.radius} fm; #it{{k}}* (MeV/#it{{c}}); SE(#it{{k}}*)'])
    hist_SE = fill_hist_from_arrays(values, errors, hist_spec, 'output/same_event.pdf')

    outfile = TFile.Open('output/correlation.root', 'RECREATE')
    hist_CF.Write()
    hist_SE.Write()
    outfile.Close()
