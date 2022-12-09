require 'histogram/array'
require 'csv'
require 'integration'
require 'distribution'
require 'pycall/import'
include PyCall::Import


module PERT
  ## 4 utility function will be used in code
  #
  # Создаем CSV file чтобы сохранить output
  puts("Make sure Results.csv file is not open")
  @csv = CSV.open("Results.csv", "w")
  def PERT.export_to_csv(bins , freqs, method_name)

    @csv << [method_name]
    @csv << bins
    @csv << freqs

  end
  def PERT.get_random(min, max)
    return rand * (max - min) + min
  end
  def PERT.pdf(x, a,b,c)
    alpha = (4*b+c-5*a)/(c-a)
    beta = (5*c-a-4*b)/(c-a)
    #puts(Math.beta(1,3))
    output = ((x-a)**(alpha-1)*(c-x)**beta-1)/(Math.beta(alpha,beta)*(c-a)**(alpha+beta+1))
  end

  def PERT.inverse_cdf(u,alpha,beta)
    # calculates inverse regularized beta function
    include PyCall::Import
    pyfrom :scipy, import: :special
    output = special.betaincinv(alpha, beta, u)
    #puts(output)
    return output
  end
  def PERT.matematics(samples)

    mean = samples.sum(0.0) / samples.size
    sum = samples.sum(0.0) { |element| (element - mean) ** 2 }
    variance = sum / (samples.size - 1)
    standard_deviation = Math.sqrt(variance)
    @csv << ["Mean", "Standard_deviation", "Variance"]
    @csv << [mean,standard_deviation,variance]
    @csv << [""]
  end




  ## Main Methods
  def PERT.inverse_method(a,b,c,number_of_experiments)
    #interval = [a,c]
    # b = mode
    step = 1.0/number_of_experiments
    alpha = (4*b+c-5*a)/(c-a)
    beta = (5*c-a-4*b)/(c-a)
    # here we ctreat a uniform ditribution in range [0,1] to sample from inverse of our cdf function
    unifrom = []
    Range.new(0,0.9999).step(step) {|x| unifrom.push(x)}
    # multiply c to all samples in order to normalize the output
    pert_samples = unifrom.map { |i| (inverse_cdf i,alpha,beta).to_f * c }
    (bins, freqs) = pert_samples.histogram(30)
    export_to_csv bins, freqs, "Inverse Method"
    puts(bins.size)
    matematics pert_samples

  end
  def PERT.metropolis_method(a,b,c,number_of_experiments)

    samples = []

    burn_in = (number_of_experiments*0.2).to_int
    # We increase number of expriments 20%, because at the end we want to remove 20% of them
    number_of_experiments = (number_of_experiments * 1.2).to_int

    # choose a random number between min and max
    max = c
    min = a
    current = PERT.get_random min,max
    for i in 1..number_of_experiments do
      samples.push(current)
      movement = PERT.get_random min,max

      curr_prob = PERT.pdf(current,a,b,c)
      move_prob = PERT.pdf(movement,a,b,c)

      acceptance = [move_prob/curr_prob,1].min
      if acceptance> rand
        current = movement
      end
    end
    # burn the initial results since they may not be so accurate
    samples = samples[burn_in..]
    (bins, freqs) = samples.histogram(30)
    PERT.export_to_csv bins, freqs,"Metropolise Method"
    matematics samples
  end
  def PERT.neyman_method(a,b,c,number_of_experiments)
    # neyman or accept and reject method
    # maximum of pdf function occurs at mode of pdf which is b here
    maximum_of_pdf = PERT.pdf(b ,a,b,c)

    samples = []
    while samples.length < number_of_experiments do
      x = PERT.get_random a, c
      y = PERT.get_random 0, maximum_of_pdf
      pdf = PERT.pdf(x,a,b,c)
      if y <= pdf
        samples.push(x)
      end


    end
    (bins, freqs) = samples.histogram(30)
    export_to_csv bins, freqs, "Neyman or Rejection Method"
    matematics samples

  end

  puts("This code samples from PERT distribution with 3 methods: Inverse, Metropolis and Neyman enter these parameters")
  puts("Enter parameters")
  print "a="
  STDOUT.flush
  a = gets.chomp.to_f
  print "b="
  STDOUT.flush
  b = gets.chomp.to_f
  while a >= b
    puts("b should be greater than a")
    print "b="
    STDOUT.flush
    b = gets.chomp.to_f
  end
  print "c="
  STDOUT.flush
  c = gets.chomp.to_f
  while c<=a or c <= b
    puts("C is the end of range. it should be greater than a and b")
    print "c="
    STDOUT.flush
    c = gets.chomp.to_f
  end
  print "Number of experiments(Samples)="
  STDOUT.flush
  n = (gets.chomp).to_i



  inverse_method a,b,c,n
  metropolis_method a,b,c,n
  neyman_method a,b,c,n

  #inverse_method 0.0,30.0,100.0,10000
  #metropolis_method 0.0,30.0,50.0,10000
  #neyman_method 0.0,30.0,50.0,10000

  puts("Results.csv file is saved")

end