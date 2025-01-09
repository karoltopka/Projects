using Microsoft.JSInterop;
using System.Text.Json;
using BlazorApp1.Shared.Enums;
using BlazorApp1.Shared.Models;
using BlazorApp1.Shared.Services;

namespace BlazorApp.client.Service
{
    internal sealed class TaskService : ITaskService
    {
        private List<TaskModel> tasks = new List<TaskModel>();
        private readonly IJSRuntime runtime;

        public TaskService(IJSRuntime runtime)
        {
            this.runtime = runtime;
        }

        public async Task Create(NewTask data)
        {
            await LoadData();
    
            var task = new TaskModel(
                Id: Guid.NewGuid(),
                Title: data.Title,
                Description: data.Description,
                Status: Status.ToDo
            );
            tasks.Add(task);

            await SaveData();
        }



        public async Task<IEnumerable<TaskModel>> GetAll(Guid? withCategory)
        {
            await LoadData();
            return tasks;
        }

        public async Task<IEnumerable<TaskModel>> GetWithStatus(Status status)
        {
            await LoadData();
            return tasks.Where(x => x.Status == status) ?? new List<TaskModel>();
        }

        public async Task ChangeStatus(Guid taskId, Status newStatus)
        {
            await LoadData();

            var task = tasks.FirstOrDefault(x => x.Id == taskId);
            if (task is null) return;

            var index = tasks.IndexOf(task);
            task = task with { Status = newStatus };
            tasks[index] = task;

            await SaveData();
        }

        public async Task Remove(Guid taskId)
        {
            await LoadData();

            var taskToRemove = tasks.FirstOrDefault(x => x.Id == taskId);
            if (taskToRemove is null) return;

            tasks.Remove(taskToRemove);
            await SaveData();
        }

        public async Task Update(UpdateTask data)
        {
            await LoadData();

            var task = tasks.FirstOrDefault(x => x.Id == data.TaskId);
            if (task is null) return;

            var index = tasks.IndexOf(task);
            task = task with { Title = data.Title, Description = data.Description };
            tasks[index] = task;

            await SaveData();
        }

        private async Task LoadData()
        {
            var loadedTasks = await runtime.InvokeAsync<string>("localStorage.getItem", "tasks");
            if (loadedTasks is not null)
            {
                tasks = JsonSerializer.Deserialize<IEnumerable<TaskModel>>(loadedTasks).ToList();
            }
        }

        private async Task SaveData()
        {
            await runtime.InvokeAsync<IEnumerable<TaskModel>>("localStorage.setItem", "tasks", JsonSerializer.Serialize(tasks));
        }

        // Wymagane przez interfejs, ale nieużywane
        public Task ChangeCategory(Guid taskId, Guid categoryId) => Task.CompletedTask;
    }
}
