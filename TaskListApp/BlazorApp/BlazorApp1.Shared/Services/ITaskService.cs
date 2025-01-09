using BlazorApp1.Shared.Enums;
using BlazorApp1.Shared.Models;

namespace BlazorApp1.Shared.Services;

public interface ITaskService
{
    Task<IEnumerable<TaskModel>> GetAll(Guid? withCategory = null);
    Task<IEnumerable<TaskModel>> GetWithStatus(Status status);
    Task Create(NewTask data);
    Task Update(UpdateTask data);
    Task Remove(Guid taskId);
    Task ChangeStatus(Guid taskId, Status newStatus);
    Task ChangeCategory(Guid taskId, Guid categoryId);
}