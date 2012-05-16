# bio-hmmer3_report

[![Build Status](https://secure.travis-ci.org/wwood/bioruby-hmmer3_report.png)](http://travis-ci.org/wwood/bioruby-hmmer3_report)

Parser class for hmmsearch and hmmscan in the HMMER 3 package.

## Examples
Input from file:

    report = Bio::Hmmer3Report.new('/path/to/hmmer_domtblout.output')
    report.hits.each do |hit|
      puts hit.target_name
      puts hit.target_accession
      puts hit.query_name
      puts hit.query_accession
      puts hit.query_length
      puts hit.full_sequence_e_value
      puts hit.full_sequence_score
      puts hit.domain_number
      puts hit.domain_sum
      puts hit.domain_c_e_value
      puts hit.domain_i_e_value
      puts hit.domain_score
      puts hit.domain_bias
      puts hit.hmm_coord_from
      puts hit.hmm_coord_to
      puts hit.ali_coord_from
      puts hit.ali_coord_to
      puts hit.env_coord_from
      puts hit.env_coord_to
      puts hit.acc
      puts hit.target_description
    end

Input from string:

    data = String.new
    data         << '                                                                           --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord'
    data << "\n" << 'target name        accession   tlen query name           accession   qlen   E-value  score  bias    of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target'
    data << "\n" << '#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------'
    data << "\n" << 'Bcl-2                PF00452.13   101 sp|P10415|BCL2_HUMAN -            239   3.7e-30  103.7   0.1   1   1   7.9e-34   4.9e-30  103.3   0.0     1   101    97   195    97   195 0.99 Apoptosis regulator proteins, Bcl-2 family'
    data << "\n" << 'BH4                  PF02180.11    27 sp|P10415|BCL2_HUMAN -            239   3.9e-15   54.6   0.1   1   1   1.3e-18   8.2e-15   53.6   0.1     2    26     8    32     7    33 0.94 Bcl-2 homology region 4'
    data << "\n"
    
    report = Bio::Hmmer3Report.new(data)
    report.hits.each do |hit|
      puts hit.target_name
      puts hit.full_sequence_e_value
    end

RDF/XML output:

    puts report.to_rdf

== References
* HMMER  http://hmmer.janelia.org/

Note: this software is under active development!

## Installation

```sh
    gem install bio-hmmer3_report
```

## Usage

```ruby
    require 'bio-hmmer3_report'
```

The API doc is online. For more code examples see the test files in
the source tree.
        
## Project home page

Information on the source tree, documentation, examples, issues and
how to contribute, see

  http://github.com/wwood/bioruby-hmmer3_report

The BioRuby community is on IRC server: irc.freenode.org, channel: #bioruby.

## Cite

If you use this software, please cite one of
  
* [BioRuby: bioinformatics software for the Ruby programming language](http://dx.doi.org/10.1093/bioinformatics/btq475)
* [Biogem: an effective tool-based approach for scaling up open source software development in bioinformatics](http://dx.doi.org/10.1093/bioinformatics/bts080)

## Biogems.info

This Biogem is published at [#bio-hmmer3_report](http://biogems.info/index.html)

## Copyright

Copyright (c) 2012 Ben J Woodcroft. See LICENSE.txt for further details.

