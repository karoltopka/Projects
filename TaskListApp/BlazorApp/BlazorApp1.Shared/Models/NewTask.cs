namespace BlazorApp1.Shared.Models;

public record NewTask(string Title, string? Description, Guid? CategoryId);