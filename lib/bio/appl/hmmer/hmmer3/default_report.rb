
require 'bio/appl/hmmer/hmmer3/tabular_report'


module Bio
  class HMMER
    class HMMER3

    def self.reports(multiple_report_text, options={})
      if [:domtblout, :tblout].include?(options[:format])
        return TabularReport.new(multiple_report_text, options[:format])
      else
        ary = []
        multiple_report_text.each_line("\n//\n") do |report|
          if block_given?
            yield DefaultHMMSearchReport.new(report)
          else
            ary << DefaultHMMSearchReport.new(report)
          end
        end
        return ary
      end
    end

    # This class is for parsing HMMSearch outputs from the default output
    class DefaultHMMSearchReport
      # Delimiter of each entry for Bio::FlatFile support.
      DELIMITER = RS = "\n//\n"

      def initialize(data)
        # The input data is divided into chunks, a hash called @report_chunks (previously called subdata)
        @report_chunks = get_subdata(data)
        
        @log = Bio::Log::LoggerPlus['bio-hmmer3report']
      end
  
      # Return the report split up into chunks, so that those chunks can be further
      # processed. Chunks are returned as a hash
      def get_subdata(data)
        subdata = {}
        header_prefix = '\Ahmmsearch :: search' ## # hmmsearch :: search profile(s) against a sequence database
        query_prefix  = '^Query:' ## Query:       2-Hacid_dh  [M=133]
        hit_prefix    = '^Scores for complete sequences' ## Scores for complete sequences (score includes all domains):
        aln_prefix    = '^Domain annotation for each sequence' ## Domain annotation for each sequence (and alignments):
        stat_prefix   = '^\nInternal pipeline statistics summary:' ## Internal pipeline statistics summary:
  
        # if header exists, get it. Header only occurs in the first report
        if data =~ /#{header_prefix}/
          subdata["header"] = data[/(\A.+?)(?=#{query_prefix})/m]
        end
  
        # split rest of report into sub-sections
        subdata["query"] =      data[/(#{query_prefix}.+?)(?=#{hit_prefix})/m]
        subdata["hit"] =        data[/(#{hit_prefix}.+?)(?=#{aln_prefix})/m]
        subdata["alignment"] =  data[/(#{aln_prefix}.+?)(?=#{stat_prefix})/m]
        subdata["statistics"] = data[/(#{stat_prefix}.+)\z/m]

        return subdata
      end
      private :get_subdata
      
      # TODO: parse header information
      # TODO: parse statistical information
      # TODO: parse sequence-wise hits
      
      # Return an array of HMMER3::Hit objects from this report
      def hits
        return [] unless @report_chunks['alignment'].match(/^>>/)
        # For each hit sequence (hits)
        sequence_annotations = @report_chunks['alignment'].split(">>")
        # puts "Found #{sequence_annotations.length} different hits e.g. #{sequence_annotations[0]}\n\n and #{sequence_annotations[1]} and \n\n #{sequence_annotations[sequence_annotations.length-1]}"
        
        alignments = []
        sequence_annotations.each_with_index do |seq_annot, i|
          #Ignore the first split as it is rubbish leftover from the split above
          next if i==0
          
          # Now split on \n\n. Each of these splits should have 1 or more domains associated
          
          #First of this split will be "stanzas" like this
          #>> 637984252  Acid345_2236 D-isomer specific 2-hydroxyacid dehydrogenase, NAD-binding [Korebacter versatilis Ellin345]
          #   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
          # ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
          #   1 !  120.2   0.5   7.6e-39   1.2e-35       1     133 []       3     313 ..       3     313 .. 0.98
          
          #And AFTER the first split comes "stanzas" like this
          #  Alignments for each domain:
          #  == domain 1    score: 120.2 bits;  conditional E-value: 7.6e-39
          #                 EEECST.-CCHHHHHCC..TEEEEEEC.GSSHHHHHC....-SEEEE-TTS-BSHHHHCC-TT--EEEES----TTB-HHHHHH---EEE--TTTTHHH CS
          #                 xxxxxxxxxxxxxxxxx..xxxxxxxx.xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx RF
          #  2-Hacid_dh   1 vlileplreeelellke..gvevevkd.ellteellekakdadalivrsntkvtaevleklpkLkviatagvGvDniDldaakerGIlVtnvpgystes 96 
          #                 +++ e++ +++++l+k+    +v++ d   ++e+lle++k+adalivrs   v+a++le++ +L+vi++agvGvDni l+aa+++GI+V+n+pg+++ +
          #   637982115   3 IVVAEKIAKAAIDLFKQdpTWNVVTPDqVAQKEQLLEQLKGADALIVRSAVFVDAAMLEHADQLRVIGRAGVGVDNIELEAATRKGIAVMNTPGANAIA 101
          #                 6889999**********8544777777778888****************************************************************** PP
          stanzas = seq_annot.split("\n\n")
          
          hit = Hit.new
          hsps = []
          
          # Parse the first
          lines = stanzas[0].split("\n")
          sequence_name = lines[0].gsub(/^ /,'')
          hit.sequence_name = sequence_name
          lines[3..(lines.length-1)].each do |line|
            annotation = Hit::Hsp.new
            
            splits = line.split(/\s+/)
            i = 1
            annotation.number = splits[i].to_i; i+=2
            annotation.score = splits[i].to_f; i+=1
            annotation.bias = splits[i].to_f; i+=1
            annotation.c_evalue = splits[i].to_f; i+=1
            annotation.i_evalue = splits[i].to_f; i+=1
            annotation.hmmfrom = splits[i].to_i; i+=1
            annotation.hmm_to = splits[i].to_i; i+=2
            annotation.alifrom = splits[i].to_i; i+=1
            annotation.ali_to = splits[i].to_i; i+=2
            annotation.envfrom = splits[i].to_i; i+=1
            annotation.env_to = splits[i].to_i; i+=2
            annotation.acc = splits[i].to_f
            hsps.push annotation
          end
          
          # Parse the second stanza and beyond
          current_hsp = nil
          (1..(stanzas.length-1)).each_with_index do |aln, index|
            next if stanzas[aln] == "\n"
            stanza = stanzas[aln].split("\n")
            line_offset = 0
            line_offset += 1 if index==0 #to account for the "Alignments for each domain:" line
            # Is this a new HSP being described?
            if matches = stanza[line_offset].match(/^  == domain (\d+)/)
              domain_index = matches[1].to_i-1
              current_hsp = hsps[domain_index]
              line_offset += 1 # to account for the "== domain 1    score: 26.8 bits;  conditional E-value: 5.7e-10" line
            end
            
            # Detect CS and RF lines
            line_offset += 1 if stanza[line_offset].split(/\s+/)[2] == 'CS'
            line_offset += 1 if stanza[line_offset].split(/\s+/)[2] == 'RF'
            
            # Add the lines to the relevant places
            current_hsp.hmmseq ||= ''
            current_hsp.hmmseq += stanza[line_offset].split(/\s+/)[3]
            current_hsp.flatseq ||= ''
            current_hsp.flatseq += stanza[2+line_offset].split(/\s+/)[3]
          end
          
          hit.hsps = hsps
          alignments.push hit
        end
        
        return alignments
      end
      
      # TODO: There is some overlapping code here between the tabular report Hit object and this object, probably should DRY it up a bit.
      class Hit
        attr_accessor :sequence_name
        
        attr_accessor :hsps
        
        def initialise
          @hsps = []
        end
        
        class Hsp
          attr_accessor :number, :score, :bias, :c_evalue, :i_evalue, :hmmfrom, :hmm_to, :alifrom, :ali_to, :envfrom, :env_to, :acc
        
          attr_accessor :hmmseq, :flatseq
        end # class DomainHitAnnotation
      end # class DomainHitAnnotation
      
      end #class DefaultHMMSearchReport
    end # class HMMER3

  end # class HMMER

end # module Bio
