import __init__
import time
from alive_progress import alive_bar

class Executor:
    
    
    def ThreadExecutor(function_or_name, *args,max_workers=''):
        import concurrent.futures
        
        # 判断函数名称是否为字符串，并尝试将其转换为函数对象
        _function = function_or_name if callable(function_or_name) else globals().get(function_or_name, None)
        if _function is None:
            raise ValueError(f"No such function: {function_or_name}")

        # 修正参数列表以适应函数需要的格式
        arguments_list = list(zip(*args)) if len(args) > 0 and len(args[0]) > 0 else [()]

        if max_workers == '':
            max_workers = min(23, len(arguments_list))
        else:
            max_workers=int(max_workers)
            
        total_tasks = len(arguments_list)
        results = [None] * total_tasks  # 初始化结果列表
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor, alive_bar(total_tasks) as bar:
            # 先提交初始批次的任务
            futures = {executor.submit(_function, *params): index for index, params in enumerate(arguments_list[:max_workers])}
            running_tasks = max_workers

            # 处理已完成的任务并提交新任务
            while futures:
                # 等待任一任务完成
                completed_future = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                for future in completed_future.done:
                    index = futures.pop(future)  # 获取完成的任务索引
                    try:
                        # 获取任务结果，并保存到结果列表中
                        results[index] = future.result()
                        bar()  # 更新进度条
                    except Exception as e:
                        print(f"Exception in task {index}: {e}")

                    # 如果还有未提交的任务，提交新任务
                    if running_tasks < total_tasks:
                        new_params = arguments_list[running_tasks]
                        new_future = executor.submit(_function, *new_params)
                        futures[new_future] = running_tasks
                        running_tasks += 1
        
        if results and isinstance(results[0], tuple):
            # 如果函数返回多个值，则将相同位置的输出合并到各自的列表中
            #num_outputs = len(results[0])
            organized_results = tuple([list(output) for output in zip(*results)])
            return organized_results
        else:
            # 如果函数返回单一值，则直接返回结果列表
            return results
       
    def ProcessExecutor(function_or_name, *args,max_workers=''):
        from concurrent.futures import ProcessPoolExecutor, as_completed
        import multiprocessing

        # 如果 function_or_name 是字符串，尝试将其转换为函数对象
        _function = function_or_name if callable(function_or_name) else globals().get(function_or_name, None)
        if _function is None:
            raise ValueError(f"找不到函数：{function_or_name}")

        # 准备参数列表
        arguments_list = list(zip(*args)) if args else [()]

        total_tasks = len(arguments_list)
        # 设置最大进程数，确保不超过任务数和处理器数
        if max_workers == '':
            max_workers = min(multiprocessing.cpu_count(), total_tasks, 4)
        else:
            max_workers=int(max_workers)


        # 初始化结果列表
        results = [None] * total_tasks

        # 使用进程池管理多进程任务
        with ProcessPoolExecutor(max_workers=max_workers) as executor, alive_bar(total_tasks) as bar:
            # 提交初始任务批次
            futures_to_index = {executor.submit(_function, *params): i for i, params in enumerate(arguments_list)}
            running_tasks = len(futures_to_index)

            # 处理完成的任务并提交新任务
            while futures_to_index:
                # 等待任一任务完成
                for future in as_completed(futures_to_index):
                    index = futures_to_index.pop(future)
                    try:
                        # 获取任务结果，并保存到结果列表
                        results[index] = future.result()
                    except Exception as e:
                        print(f"任务 {index} 中发生异常: {e}")
                    bar()  # 更新进度条

                    # 如果还有未提交的任务，提交新任务
                    if running_tasks < total_tasks:
                        new_params = arguments_list[running_tasks]
                        new_future = executor.submit(_function, *new_params)
                        futures_to_index[new_future] = running_tasks
                        running_tasks += 1

        return results
    
    def CudaExecutor(function_or_name, *args, batch_size=256, block_size=256):
        from numba import cuda
        import numpy as np
        import concurrent.futures
        import math
        """
        改进的GPU执行器，模仿ThreadExecutor的逻辑，支持动态分批调度任务到GPU。
        
        :param function: 被@cuda.jit装饰过的CUDA函数。
        :param args_list: 包含每个任务参数的列表。
        :param batch_size: 每批处理的任务数量。
        :param block_size: CUDA内核的线程块大小。
        """
        total_tasks = len(args)
        tasks_submitted = 0  # 已提交的任务计数
        results = [None] * total_tasks  # 初始化结果列表
        
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = {}

            while tasks_submitted < total_tasks or futures:
                # 计算本批次应提交的任务数
                tasks_to_submit = min(batch_size - len(futures), total_tasks - tasks_submitted)
                
                # 提交新的一批任务
                for _ in range(tasks_to_submit):
                    if tasks_submitted < total_tasks:
                        gpu_args = function_or_name(args[tasks_submitted])
                        future = executor.submit(function, *gpu_args)
                        futures[future] = tasks_submitted
                        tasks_submitted += 1
                
                # 等待任一任务完成
                done, _ = concurrent.futures.wait(futures.keys(), return_when=concurrent.futures.FIRST_COMPLETED)
                
                for future in done:
                    index = futures.pop(future)
                    try:
                        results[index] = future.result()
                    except Exception as e:
                        print(f"Exception in task {index}: {e}")
                        results[index] = None

        return results
        
    def TestProcess(function_or_name, *args):
        start_time = time.time()
        Executor.ProcessExecutor(function_or_name, *args)
        return time.time() - start_time
    
    def TestThread(function_or_name, *args):
        start_time = time.time()
        Executor.ThreadExecutor(function_or_name, *args)
        return time.time() - start_time
    
    def getmax_Thread(function, *args,max_threads):
        from Tool.Draw import Draw
        optimal_threads = 1
        shortest_time = float('inf')
        times = []

        # Test different numbers of threads
        for num_threads in range(1, max_threads + 1):
            start_time = time.time()

            # We assume ThreadExecutor has similar usage to ThreadPoolExecutor
            Executor.ThreadExecutor(function, *args, max_workers=num_threads)

            elapsed_time = time.time() - start_time
            times.append(elapsed_time)

            print(f"Time with {num_threads} threads: {elapsed_time:.2f} seconds.")

            if elapsed_time < shortest_time:
                shortest_time = elapsed_time
                optimal_threads = num_threads

        # Print times for all thread counts
        for num_threads, elapsed_time in enumerate(times, start=1):
            print(f"Threads: {num_threads}, Time: {elapsed_time:.2f} seconds")
            
        Draw.draw_line(list(range(1,max_threads+1)),times)
        return optimal_threads
        
    def getmax_Process(function, *args,max_threads):
        from Tool.Draw import Draw
        optimal_threads = 1
        shortest_time = float('inf')
        times = []

        # Test different numbers of threads
        for num_threads in range(1, max_threads + 1):
            start_time = time.time()

            # We assume ThreadExecutor has similar usage to ThreadPoolExecutor
            Executor.ProcessExecutor(function, *args, max_workers=num_threads)

            elapsed_time = time.time() - start_time
            times.append(elapsed_time)

            print(f"Time with {num_threads} threads: {elapsed_time:.2f} seconds.")

            if elapsed_time < shortest_time:
                shortest_time = elapsed_time
                optimal_threads = num_threads

        # Print times for all thread counts
        for num_threads, elapsed_time in enumerate(times, start=1):
            print(f"Threads: {num_threads}, Time: {elapsed_time:.2f} seconds")
            
        Draw.draw_line(list(range(1,max_threads+1)),times)
        return optimal_threads
              
    def chooseExecutor(function_or_name, *args, sample_size=100):
        # 提取样本任务
        sample_args = [arg[:sample_size] for arg in args]

        # 测试多进程执行器
        process_time=Executor.TestProcess(function_or_name, *sample_args)

        # 测试多线程执行器
        thread_time = Executor.TestThread(function_or_name, *sample_args)

        # 比较执行时间并选择
        if process_time < thread_time:
            print(f"多进程更快，耗时 {process_time} 秒，将使用多进程执行器。")
            return Executor.ProcessExecutor
        else:
            print(f"多线程更快，耗时 {thread_time} 秒，将使用多线程执行器。")
            return Executor.ThreadExecutor


    def GaussianTgread(function_or_name, *args,max_workers=4, total_cores=8):
        import concurrent.futures
        # 转换函数名称为函数对象
        _function = function_or_name if callable(function_or_name) else globals().get(function_or_name, None)
        if _function is None:
            raise ValueError(f"No such function: {function_or_name}")

        # 准备参数列表
        arguments_list = list(zip(*args))
        total_tasks = len(arguments_list)
        
        # 初始化结果列表
        results = [None] * total_tasks

        # 初始化任务管理
        task_status = [{'core_id': None, 'is_running': False} for _ in range(total_tasks)]
        available_cores = list(range(total_cores))

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            with alive_bar(total_tasks) as bar:
                futures_to_task = {}

                # 提交初始批次的任务
                for i in range(min(max_workers, total_tasks)):
                    if available_cores:
                        core_id = available_cores.pop(0)
                        task_status[i]['core_id'] = core_id
                        task_status[i]['is_running'] = True
                        future = executor.submit(_function, core_id, *arguments_list[i])
                        futures_to_task[future] = i

                # 动态管理任务
                while futures_to_task:
                    # 等待任一任务完成
                    done, _ = concurrent.futures.wait(futures_to_task, return_when=concurrent.futures.FIRST_COMPLETED)

                    for future in done:
                        task_index = futures_to_task.pop(future)
                        task_status[task_index]['is_running'] = False
                        core_id = task_status[task_index]['core_id']
                        available_cores.append(core_id)  # 将核心ID放回可用列表

                        try:
                            # 获取任务结果，并保存到结果列表中
                            results[task_index] = future.result()
                            bar()  # 更新进度条
                        except Exception as e:
                            print(f"Exception in task {task_index}: {e}")

                        # 如果还有未提交的任务，则提交新任务
                        for next_task_index, status in enumerate(task_status):
                            if not status['is_running'] and available_cores and results[next_task_index] is None:
                                next_core_id = available_cores.pop(0)
                                task_status[next_task_index]['core_id'] = next_core_id
                                task_status[next_task_index]['is_running'] = True
                                next_future = executor.submit(_function, next_core_id, *arguments_list[next_task_index])
                                futures_to_task[next_future] = next_task_index
                                break  # 一次只提交一个新任务，以保持最大并发数


        return results