# HPC - `CUDA project`

## Installation

Run the following commands to install the project:
```bash
make
```

In order to run the project you can use:
```bash
make run
```

To generate the report you can use:
```bash
make report
```

To clean the project you can use:
```bash
make clean
```

## `vep` install process

```bash
sudo cpan Archive::Zip DBI Module::Build
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl
```

