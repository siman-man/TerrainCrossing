class Analyze
  TSV_LABEL = ['seed', 'score', 'time'].join("\t")

  def initialize
    data_list = parse
    time_stamp = Time.now.strftime("%m%d-%H%m%S")

    @record = File.open(File.expand_path("data/record-#{time_stamp}.tsv", Dir::pwd), 'w')

    begin
      @record.puts(TSV_LABEL)

      sum_score = 0.0
      seeds = []

      data_list.each do |data|
        data = clean_data(data)

        seeds << data[0]
        puts data[1]
        sum_score += data[1]
        @record.puts(data.join("\t"))
      end

      #p seeds
      puts sum_score
    ensure
      @record&.close
    end
  end 

  def clean_data(data)
    seed = data['seed'].to_i
    score = data['score'].to_f
    /(?<minute>\d+)m(?<second>(\d|\.)+)s/ =~ data['user']
    time = minute.to_f * 60 + second.to_f

    [seed, score, time]
  end

  def parse
    filepath = File.expand_path('result.txt', Dir::pwd)
    data_list = []
    data = {}

    File.open(filepath, 'r') do |file|
      check = false
      file.each_line do |line|
        line = clean_line(line)

        if line =~ /begin/
          data = {}
          check = false
        elsif line =~ /!end!/
        #elsif line =~ /!end!/ && check
          data_list << data.dup
        else
          if validate(line)
            h = line2hash(line)
            data.merge!(h)
          end
          if validate2(line)
            #puts line
            check ||= true
          end
        end
      end
    end

    data_list
  end

  def line2hash(line)
    Hash[*line.split]
  end

  def clean_line(line)
    line.chomp.downcase.delete('=')
  end

  def validate(line)
    line =~ /^(score|seed|user)/
  end

  def validate2(line)
    line =~ /nship  3 nufo  [1-9]/
  end
end

Analyze.new
