using Downloads

const ROOT = joinpath(@__DIR__, "..")

# Pre-selected VLBI clean model FITS images from the Astrogeo Center (NASA mirror).
# One file per source: highest-frequency band available, Gaia DR3 era preferred,
# observation with most clean components.
const FILENAMES = [
    "J0045+2127_X_2010_03_23_pet_map.fits",
    "J0113+4948_X_2016_07_22_pet_map.fits",
    "J0222+4302_U_2003_06_21_moj_map.fits",
    "J0225-2312_X_2017_09_18_pet_map.fits",
    "J0303-2407_U_2012_12_10_moj_map.fits",
    "J0521+2112_U_2014_01_11_moj_map.fits",
    "J0721+7120_U_2011_05_21_moj_map.fits",
    "J0738+1742_U_2008_11_26_moj_map.fits",
    "J0746-1555_X_2017_08_05_pet_map.fits",
    "J0757+0956_U_2025_01_06_moj_map.fits",
    "J0856-1105_U_2015_12_19_moj_map.fits",
    "J0921+6215_U_2023_08_18_moj_map.fits",
    "J1058+0133_U_2002_09_20_moj_map.fits",
    "J1221+2813_U_2009_05_14_moj_map.fits",
    "J1459+7140_U_2025_01_26_moj_map.fits",
    "J1543+0452_U_2017_03_11_moj_map.fits",
    "J1716+6836_U_2016_09_17_moj_map.fits",
    "J1728+5013_U_2015_12_19_moj_map.fits",
    "J1734+3857_U_2010_11_20_moj_map.fits",
    "J1800+7828_U_2010_12_18_moj_map.fits",
    "J1801+4404_U_2009_12_17_moj_map.fits",
    "J1806+6949_U_2019_11_14_moj_map.fits",
    "J1824+5651_U_2002_08_02_moj_map.fits",
    "J1927+6117_U_2015_10_02_moj_map.fits",
    "J2243-2544_X_2017_08_05_pet_map.fits",
]

outdir = joinpath(ROOT, "data/vlbi_images")
mkpath(outdir)

for fname in FILENAMES
    outpath = joinpath(outdir, fname)
    if isfile(outpath)
        println("  $fname (exists)")
        continue
    end
    j2000 = first(split(fname, '_'))
    url = "https://astrogeo.smce.nasa.gov/images/$j2000/$fname"
    print("  $fname ... ")
    Downloads.download(url, outpath)
    println("ok")
end

println("Downloaded $(length(FILENAMES)) FITS files to $outdir")
