#
# = bio/appl/hmmer/hmmer3/report.rb - hmmscan/hmmsearch parser
#
# Copyright::   Copyright (C) 2011
#               Christian Zmasek <cmzmasek@yahoo.com>
# License::     The Ruby License
#
# $Id:$
#

require 'stringio'

module Bio
  class HMMER
    # == Description
    #
    # Parser class for hmmsearch and hmmscan in the HMMER 3 package.
    #
    # == Examples
    #
    #    # Input from file:
    #    report = Bio::HMMER::HMMER3::Report.new('/path/to/hmmer_domtblout.output')
    #    report.hits.each do |hit|
    #      puts hit.target_name
    #      puts hit.target_accession
    #      puts hit.query_name
    #      puts hit.query_accession
    #      puts hit.query_length
    #      puts hit.full_sequence_e_value
    #      puts hit.full_sequence_score
    #      puts hit.domain_number
    #      puts hit.domain_sum
    #      puts hit.domain_c_e_value
    #      puts hit.domain_i_e_value
    #      puts hit.domain_score
    #      puts hit.domain_bias
    #      puts hit.hmm_coord_from
    #      puts hit.hmm_coord_to
    #      puts hit.ali_coord_from
    #      puts hit.ali_coord_to
    #      puts hit.env_coord_from
    #      puts hit.env_coord_to
    #      puts hit.acc
    #      puts hit.target_description
    #    end
    #
    #
    #
    #    # Input from string:
    #    data = String.new
    #    data         << '#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord'
    #    data << "\n" << '# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target'
    #    data << "\n" << '#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------'
    #    data << "\n" << 'Bcl-2                PF00452.13   101 sp|P10415|BCL2_HUMAN -            239   3.7e-30  103.7   0.1   1   1   7.9e-34   4.9e-30  103.3   0.0     1   101    97   195    97   195 0.99 Apoptosis regulator proteins, Bcl-2 family'
    #    data << "\n" << 'BH4                  PF02180.11    27 sp|P10415|BCL2_HUMAN -            239   3.9e-15   54.6   0.1   1   1   1.3e-18   8.2e-15   53.6   0.1     2    26     8    32     7    33 0.94 Bcl-2 homology region 4'
    #    data << "\n"
    #
    #    report = Bio::HMMER::HMMER3::Report.new(data)
    #    report.hits.each do |hit|
    #      puts hit.target_name
    #      puts hit.full_sequence_e_value
    #    end
    #
    #
    #
    # == References
    #
    # * HMMER
    #   http://hmmer.janelia.org/
    class HMMER3
      class Report
        def initialize(hmmer_output, format = nil)
          
          @hits = Array.new
          @line_number = 0
          @format = format
          if File.exists?(hmmer_output.to_s)
            my_hmmer_output_file = File.new(hmmer_output.to_s)
            my_hmmer_output_file.each_line() { |line| parse_line(line) }
          else
            str = StringIO.new(hmmer_output)
            str.each_line() { |line| parse_line(line) }
          end
        end
    
        attr_reader :hits
        attr_reader :format
    
        private
    
        def parse_line(line)
          @line_number += 1
          if line  =~ /^#.+this\s+domain/
            @format = :domtblout
          elsif line =~ /^#.+best\s+1\s+domain/
            @format = :tblout
          elsif line =~ /\S/ && line !~ /^#/
            if @format == nil
              if looks_like_per_domain_result?(line)
                @format = :domtblout
              else
                @format = :tblout
              end
            end
            if @format == :domtblout
              @hits << PerDomainHit.new(line, @line_number)
            elsif @format == :tblout
              @hits << PerSequenceHit.new(line, @line_number)
            else
              raise ArgumentError, "attempt to parse hmmscan/hmmsearch output style other than \"domtblout\" or \"tblout\""
            end
          end
        end
    
        def looks_like_per_domain_result?(line)
          line =~ /^(\S*)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s*(.*)/
        end
    
      end # class Report
    
      class Hit
        def initialize
          # This is an abstract class. Prevents 'new' being called on this class
          # and force implementation of 'initialize' in inheriting classes.
          raise NotImplementedError
        end
        attr_reader :target_name
        attr_reader :target_accession
        attr_reader :target_description
        attr_reader :query_name
        attr_reader :query_accession
        attr_reader :full_sequence_e_value
        attr_reader :full_sequence_score
        attr_reader :full_sequence_bias
    
      end # class Hit
    
    
      class PerSequenceHit < Hit
    
        # Sets hit data.
        def initialize(line, line_number)
    
          # tblout:
          #               tn    tacc      qn    qacc fs_eval fs_scor fs_bias   bst_e bst_scor bst_bias   exp     reg     clu      ov     env      dom     rep     inc   desc
          #                1       2       3       4       5       6       7       8       9      10      11      12      13      14      15       16      17      18     19
          if  line =~ /^(\S*)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(.*)/
            @target_name = $1
            @target_accession = $2
            @query_name = $3
            @query_accession = $4
            @full_sequence_e_value = $5.to_f
            @full_sequence_score = $6.to_f
            @full_sequence_bias =  $7.to_f
            @best_1_domain_e_value = $8.to_f
            @best_1_domain_score = $9.to_f
            @best_1_domain_bias =  $10.to_f
            @domain_number_est_exp = $11.to_i
            @domain_number_est_reg = $12.to_i
            @domain_number_est_clu = $13.to_i
            @domain_number_est_ov  = $14.to_i
            @domain_number_est_env = $15.to_i
            @domain_number_est_dom = $16.to_i
            @domain_number_est_rep = $17.to_i
            @domain_number_est_inc = $18.to_i
            @target_description = $19
          else
            raise ArgumentError, "line "+ line_number.to_s + " is in an unrecognized format [#{line}]"
          end
    
        end # initialize
    
        attr_reader :best_1_domain_e_value
        attr_reader :best_1_domain_score
        attr_reader :best_1_domain_bias
        attr_reader :domain_number_est_exp
        attr_reader :domain_number_est_reg
        attr_reader :domain_number_est_clu
        attr_reader :domain_number_est_ov
        attr_reader :domain_number_est_env
        attr_reader :domain_number_est_dom
        attr_reader :domain_number_est_rep
        attr_reader :domain_number_est_inc
    
      end # class PerSequenceHit
    
      class PerDomainHit < Hit
    
        # Sets hit data.
        def initialize(line, line_number)
    
          # domtblout:
          #                tn     acc    tlen   query     acc    qlen  Evalue   score    bias      #       of     c-E    i-E     score   bias      hf      ht      af      at     ef      et     acc  desc
          #                 1       2       3       4       5       6       7       8       9      10      11      12     13      14      15       16      17      18      19     20      21      22     23
          if  line =~ /^(\S*)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s*(.*)/
            @target_name = $1
            @target_accession = $2
            @target_length = $3.to_i
            @query_name = $4
            @query_accession = $5
            @query_length = $6.to_i
            @full_sequence_e_value = $7.to_f
            @full_sequence_score = $8.to_f
            @full_sequence_bias =  $9.to_f
            @domain_number = $10.to_i
            @domain_sum = $11.to_i
            @domain_c_e_value = $12.to_f
            @domain_i_e_value = $13.to_f
            @domain_score = $14.to_f
            @domain_bias = $15.to_f
            @hmm_coord_from = $16.to_i
            @hmm_coord_to = $17.to_i
            @ali_coord_from = $18.to_i
            @ali_coord_to = $19.to_i
            @env_coord_from = $20.to_i
            @env_coord_to = $21.to_i
            @acc = $22.to_f
            @target_description = $23
          else
            raise ArgumentError, "line "+ line_number.to_s + " is in a unrecognized format [#{line}]"
          end
    
        end # initialize
    
        attr_reader :target_length
        attr_reader :query_length
        attr_reader :domain_number
        attr_reader :domain_sum
        attr_reader :domain_c_e_value
        attr_reader :domain_i_e_value
        attr_reader :domain_score
        attr_reader :domain_bias
        attr_reader :hmm_coord_from
        attr_reader :hmm_coord_to
        attr_reader :ali_coord_from
        attr_reader :ali_coord_to
        attr_reader :env_coord_from
        attr_reader :env_coord_to
        attr_reader :acc
    
      end # class PerDomainHit
    end # class HMMER3
  end # class HMMER
end # module Bio