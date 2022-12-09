require 'histogram/array'
require 'csv'
require 'integration'
require 'distribution'
require 'pycall/import'
require 'benchmark'
include PyCall::Import


module PERT

  # Создаем CSV file чтобы сохранить output
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
    # вычисляем обратную бета-функцию
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




  ## Основные методы
  def PERT.inverse_method(a,b,c,number_of_experiments)
    #интервал = [a,c]
    step = 1.0/number_of_experiments
    alpha = (4*b+c-5*a)/(c-a)
    beta = (5*c-a-4*b)/(c-a)
    # здесь мы обрабатываем равномерное распределение в диапазоне [0,1] для выборки из обратной функции cdf
    unifrom = []
    Range.new(0,0.9999).step(step) {|x| unifrom.push(x)}
    pert_samples = unifrom.map { |i| (inverse_cdf i,alpha,beta).to_f * c }
    (bins, freqs) = pert_samples.histogram(30)
    export_to_csv bins, freqs, "Inverse Method"
    puts(bins.size)
    matematics pert_samples

  end
  def PERT.metropolis_method(a,b,c,number_of_experiments)

    samples = []

    burn_in = (number_of_experiments*0.2).to_int
    #  Увеличиваем количество экспериментов на 20%, потому что в конце хотим убрать 20% из них
    number_of_experiments = (number_of_experiments * 1.2).to_int

    # выбираем рандомное число между мин и макс
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
    # убираем первоначальные результаты, так как они могут быть не точными
    samples = samples[burn_in..]
    (bins, freqs) = samples.histogram(30)
    PERT.export_to_csv bins, freqs,"Metropolise Method"
    matematics samples
  end
  def PERT.neyman_method(a,b,c,number_of_experiments)
    # метод Неймана
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

  puts ("\nІнформація щодо використання програмного забезпечення, яке генерує випадкові велечини за PERT розподілом.")
  puts "\nУ даній роботі описано методи та алгоритми генерування випадкової величині, такі як Метод зворотної функції, Метод Неймана та Метод Метрополіса."
  puts "\nЩоб розпочати, необхідно ввести наступні значення:"
  puts ("
a - мінімальне значення якого може набувати змінна
b - найбільш ймовірне значення якого може набувати змінна (Потрібно вводити число більше, ніж а)
c - максимальне значення якого може набувати змінна (Тож потрібно вводити число більше, ніж а та b)
n - кількість експериментів")
  puts("Після коректного введення усіх даних будуть розраховані випадкові величини, і також будуть визначені час виконання,
математичне очікування, дисперсія і середньоквадратичне відхилення для кожного методу.")
  puts "\nРезультати обчислень можна буде побачити у файлі Results.csv та гістограми у файлі Output.xlsx"
  puts("\nРоботу виконала студентка групи КС-43 Алексеєнкова Дарья\n")
  puts ("\n")
  puts("Введіть параметри:")
  print "a ="
  STDOUT.flush
  a = gets.chomp.to_f
  print "b ="
  STDOUT.flush
  b = gets.chomp.to_f
  while a >= b
    puts("Параметр b повинен бути більше, ніж a")
    print "b ="
    STDOUT.flush
    b = gets.chomp.to_f
  end
  print "c ="
  STDOUT.flush
  c = gets.chomp.to_f
  while c<=a or c <= b
    puts("Параметр с повинен бути більше , ніж  a та b")
    print "c ="
    STDOUT.flush
    c = gets.chomp.to_f
  end
  print "Кількість експериментів ="
  STDOUT.flush
  n = (gets.chomp).to_i

  #inverse_method a,b,c,n
  #metropolis_method a,b,c,n
  #neyman_method a,b,c,n

  #inverse_method 0.0,30.0,100.0,10000
  #metropolis_method 0.0,30.0,50.0,10000
  #neyman_method 0.0,30.0,50.0,10000
  time = Benchmark.measure {
    puts("Відповідь:")
    metropolis_method a, b, c, n
    neyman_method a, b, c, n
    inverse_method a, b, c, n
  }
  puts("Час виконання методів (с):",time.real)
  cmd = "Output.xlsx"
  system('start "" ' + cmd)
  puts("Results.csv file is saved")

end